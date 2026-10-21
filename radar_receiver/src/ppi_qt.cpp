/*
 * ppi_qt.cpp - interactive Qt PPI (plan position indicator) of the radar receiver.
 *
 * The worker thread posts antenna positions and plots with ppi_post_antenna()
 * and ppi_post_plot(). These calls only store the data under a mutex; a timer
 * in the GUI thread (25 Hz) takes the new data and repaints the widget.
 *
 * Interaction:
 *   mouse wheel            zoom in/out around the cursor (x1 ... x500)
 *   left drag              pan
 *   double click           reset the view (whole instrumented range)
 *   hover over a plot      data tip: range, azimuth, elevation, radial velocity,
 *                          SNR, time, sub-pulse, hits, beam/cell/bin, pulse number
 *   click on a plot        pin its data tip on the display (click again: remove)
 *   right click            options: short/long-pulse plots, range rings, plot
 *                          labels, antenna trail, which plots are shown (until
 *                          the antenna passes again / 2 scans / all), minimum
 *                          SNR, clear plots, remove data tips, reset view,
 *                          pause, save as PNG
 *   keys                   + - zoom, arrows pan, R/Home reset, C remove the
 *                          data tips, Space pause/resume, S save as PNG
 * The status line shows the range, azimuth and east/north of the cursor.
 *
 * Plot persistence (default): a plot stays until the antenna sweeps over its
 * azimuth again in a later scan, like on a real PPI. For this the antenna
 * rotation is accumulated ("sweep", degrees turned since the first update), so
 * the scan boundary (the start azimuth of the antenna) does not matter.
 *
 * The file does not use Q_OBJECT, so it builds with or without moc, with
 * Qt 5 (>= 5.12) and Qt 6.
 */
#include "radar_rx.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QColor>
#include <QConicalGradient>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QFont>
#include <QFontMetrics>
#include <QKeyEvent>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QPixmap>
#include <QPointF>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QTimer>
#include <QToolTip>
#include <QWheelEvent>
#include <QWidget>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <deque>
#include <functional>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

namespace {

constexpr double kPi          = 3.14159265358979323846;
constexpr double kDeg         = kPi / 180.0;
constexpr size_t kMaxPlots    = 4000;     /* plots kept on the display              */
constexpr size_t kMaxPending  = 20000;    /* plots waiting for the GUI thread       */
constexpr int    kRefreshMs   = 40;       /* GUI refresh period (25 Hz)             */
constexpr double kHeaderPx    = 48.0;     /* status text above the display          */
constexpr double kLabelPx     = 28.0;     /* room for the azimuth labels            */
constexpr double kPlotRadius  = 3.5;      /* plot symbol radius [px]                */
constexpr double kPickPx      = 8.0;      /* data tip: max. distance to a plot [px] */
constexpr double kTrailDeg    = 20.0;     /* fading trail behind the antenna line   */
constexpr double kMaxZoom     = 500.0;

const QColor kColorShort(0, 230, 255);    /* cyan                                   */
const QColor kColorLong(255, 230, 0);     /* yellow                                 */
const QColor kColorAntenna(0, 220, 0);    /* green                                  */
const QColor kColorGrid(70, 70, 70);
const QColor kColorRim(120, 120, 120);
const QColor kColorLabel(150, 150, 150);
const QColor kColorPin(80, 230, 80);

double wrap360(double a)
{
    a = std::fmod(a, 360.0);
    return a < 0.0 ? a + 360.0 : a;
}

double wrap180(double a)
{
    a = wrap360(a);
    return a >= 180.0 ? a - 360.0 : a;
}

/* ------------------------------------------------------------------------ */
/* Data shared between the worker thread and the GUI thread                  */
/* ------------------------------------------------------------------------ */

struct AntennaState {
    bool     valid = false;
    double   az_deg = 0.0;
    double   time_s = 0.0;
    uint32_t scan = 0;
    double   sweep_deg = 0.0;   /* rotation accumulated since the first update */
};

struct PlotMark {
    double    east_m = 0.0;     /* ground position                              */
    double    north_m = 0.0;
    int       pulse = PULSE_SHORT;
    double    sweep_deg = 0.0;  /* antenna sweep at which it pointed at the plot */
    uint32_t  scan = 0;         /* antenna scan when the plot arrived           */
    uint64_t  id = 0;           /* unique number (pinned data tips)             */
    Detection det;              /* the plot itself (data tip)                   */
};

struct Shared {
    std::mutex            mutex;
    AntennaState          antenna;
    std::vector<PlotMark> newPlots;   /* plots not yet taken by the GUI thread */
    bool                  clear = false;   /* a new run started: clear the display */
    uint64_t              nextId = 1;
};

Shared &shared()
{
    static Shared s;
    return s;
}

void reset_shared()
{
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.antenna = AntennaState();
    s.newPlots.clear();
    s.clear = true;
}

/* Rotation from the previous antenna update to az_deg. The antenna turns
 * clockwise, so a gap in the data (e.g. sector transmission) is a step forward;
 * a tiny step backwards (jitter) stays a step backwards. When the scan number
 * jumped by more than one, the missing whole turns are added. */
double antenna_step(const AntennaState &prev, double az_deg, uint32_t scan)
{
    double step = wrap360(az_deg - prev.az_deg);
    if (step > 350.0) step -= 360.0;
    const long long jump = (long long)scan - (long long)prev.scan;
    if (jump > 1) step += 360.0 * (double)(jump - 1);
    return step;
}

/* Ring spacing for the visible range: about 3 to 6 rings. */
double ring_step_m(double range_m)
{
    const double km = range_m / 1e3;
    if (km > 500.0) return 100e3;
    if (km > 200.0) return 50e3;
    if (km > 60.0)  return 20e3;
    if (km > 30.0)  return 10e3;
    if (km > 15.0)  return 5e3;
    if (km > 6.0)   return 2e3;
    if (km > 3.0)   return 1e3;
    if (km > 1.5)   return 500.0;
    if (km > 0.6)   return 200.0;
    return 100.0;
}

QString km_text(double m)
{
    return QString::number(m / 1e3, 'g', 4) + QStringLiteral(" km");
}

QString deg_text(double deg, int decimals)
{
    return QString::number(deg, 'f', decimals) + QChar(0x00B0);
}

QFont pixel_font(const QFont &base, int pixels, bool bold = false)
{
    QFont f(base);
    f.setPixelSize(pixels);
    f.setBold(bold);
    return f;
}

/* Data tip text of a plot. */
QString plot_text(const Detection &d)
{
    const QString el = std::isfinite(d.el_deg) ? deg_text(d.el_deg, 2) : QStringLiteral("- (no DBF)");
    QStringList rows;
    rows << QStringLiteral("<b>Plot, %1 pulse</b>").arg(d.pulse == PULSE_LONG ? "long" : "short")
         << QStringLiteral("Range: %1 km").arg(d.range_m / 1e3, 0, 'f', 3)
         << QStringLiteral("Azimuth: %1").arg(deg_text(d.az_deg, 2))
         << QStringLiteral("Elevation: %1").arg(el)
         << QStringLiteral("Radial velocity: %1 m/s").arg(d.vr_mps, 0, 'f', 2)
         << QStringLiteral("SNR: %1 dB").arg(d.snr_dB, 0, 'f', 1)
         << QStringLiteral("Time: %1 s").arg(d.time_s, 0, 'f', 3)
         << QStringLiteral("CPI hits: %1").arg(d.hits)
         << QStringLiteral("Beam %1, range cell %2, Doppler bin %3").arg(d.beam).arg(d.rIdx).arg(d.dIdx)
         << QStringLiteral("Pulse number %1").arg(d.pulseNumber);
    return rows.join(QStringLiteral("<br>"));
}

/* Event positions (the accessors differ between Qt 5 and Qt 6). */
QPointF event_pos(const QMouseEvent *e)
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    return e->position();
#else
    return e->localPos();
#endif
}

QPointF event_pos(const QWheelEvent *e)
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    return e->position();
#else
    return e->posF();
#endif
}

/* ------------------------------------------------------------------------ */
/* PPI widget                                                                */
/* ------------------------------------------------------------------------ */

enum class History { UntilSwept, TwoScans, All };

class PpiWidget : public QWidget {
public:
    PpiWidget(double maxRange_m, volatile sig_atomic_t *stopFlag)
        : maxRange_m_(maxRange_m > 0.0 ? maxRange_m : 60e3), stopFlag_(stopFlag)
    {
        setAttribute(Qt::WA_OpaquePaintEvent);
        setMouseTracking(true);                 /* hover data tips without a button */
        setFocusPolicy(Qt::StrongFocus);        /* keys */
        setMinimumSize(300, 300);
    }

    void refresh();                     /* takes the posted data (GUI thread) */
    void setFinished(int workerCode)
    {
        finished_ = true;
        workerCode_ = workerCode;
        update();
    }

    /* view and options, also used by the tests */
    void   zoomAt(double factor, QPointF screen);
    void   panPixels(QPointF delta);
    void   resetView();
    double zoom() const { return zoom_; }
    QPointF viewCentre_m() const { return QPointF(viewE_, viewN_); }
    QPointF screenOf(double east_m, double north_m) { layoutDisplay(); return toScreen(east_m, north_m); }
    int    pick(QPointF screen) const;  /* index of the plot under the cursor, -1 if none */
    QString tipText(QPointF screen) const
    {
        const int i = pick(screen);
        return i < 0 ? QString() : plot_text(plots_[(size_t)i].det);
    }
    void   clickAt(QPointF screen);     /* pin / unpin the data tip of a plot */
    size_t pinnedCount() const { return pinned_.size(); }
    size_t plotCount() const { return plots_.size(); }
    size_t visibleCount() const;
    QString cursorText() const { return cursorText_; }
    void   setShowShort(bool on) { showShort_ = on; update(); }
    void   setShowLong(bool on) { showLong_ = on; update(); }
    void   setMinSnr(double dB) { minSnr_dB_ = dB; update(); }
    void   setHistory(History h) { history_ = h; refresh(); }
    void   setPaused(bool on) { paused_ = on; update(); }
    bool   paused() const { return paused_; }
    void   clearPlots() { plots_.clear(); pinned_.clear(); update(); }
    bool   saveSnapshot(const QString &path);

protected:
    void paintEvent(QPaintEvent *event) override;
    void closeEvent(QCloseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void leaveEvent(QEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;

private:
    void   layoutDisplay();
    void   drawBackground();
    void   drawAntenna(QPainter &p) const;
    void   drawPlots(QPainter &p) const;
    void   drawPins(QPainter &p) const;
    void   drawStatus(QPainter &p) const;
    void   drawLegend(QPainter &p) const;
    bool   expired(const PlotMark &m) const;
    bool   shown(const PlotMark &m) const;
    void   updateCursor(QPointF screen);
    void   setView(double zoom, double e, double n);
    double turned(const PlotMark &m) const { return std::fabs(antenna_.sweep_deg - m.sweep_deg); }
    double scale() const { return radius_px_ * zoom_ / maxRange_m_; }   /* pixels per metre */
    QPointF toScreen(double east_m, double north_m) const
    {
        const double k = scale();
        return QPointF(centre_.x() + k * (east_m - viewE_), centre_.y() - k * (north_m - viewN_));
    }
    QPointF toGround(QPointF s) const
    {
        const double k = scale();
        return QPointF(viewE_ + (s.x() - centre_.x()) / k, viewN_ - (s.y() - centre_.y()) / k);
    }
    QPointF polar(double az_deg, double r_m) const
    {
        return toScreen(r_m * std::sin(az_deg * kDeg), r_m * std::cos(az_deg * kDeg));
    }

    double                maxRange_m_;
    volatile sig_atomic_t *stopFlag_;
    AntennaState          antenna_;
    std::deque<PlotMark>  plots_;
    std::vector<PlotMark> incoming_;    /* buffer swapped with Shared::newPlots */
    bool                  finished_ = false;
    int                   workerCode_ = 0;
    QPointF               centre_;      /* screen point of the view centre       */
    double                radius_px_ = 1.0;   /* screen radius of maxRange at zoom 1 */
    QRectF                viewRect_;    /* drawing area below the status lines   */
    QPixmap               background_;  /* rings, spokes and labels            */
    bool                  backgroundDirty_ = true;
    QSize                 backgroundSize_;

    /* view */
    double                zoom_ = 1.0;
    double                viewE_ = 0.0, viewN_ = 0.0;   /* ground point at the view centre */
    /* mouse */
    bool                  dragging_ = false, dragMoved_ = false;
    QPointF               dragStart_;
    double                dragE_ = 0.0, dragN_ = 0.0;
    uint64_t              hoverId_ = 0;
    QString               cursorText_;
    std::set<uint64_t>    pinned_;
    /* options */
    bool                  showShort_ = true, showLong_ = true, showRings_ = true;
    bool                  showLabels_ = false, showTrail_ = true, paused_ = false;
    double                minSnr_dB_ = -1e9;
    History               history_ = History::UntilSwept;
};

void PpiWidget::refresh()
{
    if (paused_) return;                /* the data waits in Shared */
    bool clear = false;
    incoming_.clear();
    {
        Shared &s = shared();
        std::lock_guard<std::mutex> lock(s.mutex);
        incoming_.swap(s.newPlots);     /* the worker keeps the old capacity */
        antenna_ = s.antenna;
        clear = s.clear;
        s.clear = false;
    }
    if (clear) { plots_.clear(); pinned_.clear(); }
    plots_.insert(plots_.end(), incoming_.begin(), incoming_.end());
    plots_.erase(std::remove_if(plots_.begin(), plots_.end(),
                                [this](const PlotMark &m) { return expired(m); }),
                 plots_.end());
    while (plots_.size() > kMaxPlots) plots_.pop_front();
    for (auto it = pinned_.begin(); it != pinned_.end();) {   /* tips of removed plots */
        const uint64_t id = *it;
        const bool alive = std::any_of(plots_.begin(), plots_.end(), [id](const PlotMark &m) { return m.id == id; });
        it = alive ? std::next(it) : pinned_.erase(it);
    }
    update();
}

/* Default: a plot disappears when the antenna has turned a full revolution
 * since it pointed at the plot, i.e. it sweeps over the plot azimuth in the
 * next scan (safety net for gaps in the data: two scans later it is gone in
 * any case). Options: two revolutions, or all plots until cleared. */
bool PpiWidget::expired(const PlotMark &m) const
{
    if (history_ == History::All) return false;
    const double turns = history_ == History::TwoScans ? 2.0 : 1.0;
    const long long scans = (long long)antenna_.scan - (long long)m.scan;
    return turned(m) >= 360.0 * turns || scans >= 1 + (long long)turns || scans <= -2;
}

bool PpiWidget::shown(const PlotMark &m) const
{
    return (m.pulse == PULSE_LONG ? showLong_ : showShort_) && !(m.det.snr_dB < minSnr_dB_);
}

size_t PpiWidget::visibleCount() const
{
    return (size_t)std::count_if(plots_.begin(), plots_.end(), [this](const PlotMark &m) { return shown(m); });
}

int PpiWidget::pick(QPointF screen) const
{
    int best = -1;
    double bestD = kPickPx;
    for (size_t i = 0; i < plots_.size(); ++i) {
        const PlotMark &m = plots_[i];
        if (!shown(m)) continue;
        const QPointF d = toScreen(m.east_m, m.north_m) - screen;
        const double dist = std::hypot(d.x(), d.y());
        if (dist <= bestD) { bestD = dist; best = (int)i; }
    }
    return best;
}

void PpiWidget::clickAt(QPointF screen)
{
    const int i = pick(screen);
    if (i < 0) return;
    const uint64_t id = plots_[(size_t)i].id;
    if (!pinned_.erase(id)) pinned_.insert(id);
    update();
}

bool PpiWidget::saveSnapshot(const QString &path)
{
    const QPixmap pixmap = grab();
    if (pixmap.save(path, "PNG")) return true;
    std::fprintf(stderr, "PPI: cannot save the snapshot to %s\n", path.toLocal8Bit().constData());
    return false;
}

/* ---- view ---------------------------------------------------------------- */

void PpiWidget::setView(double zoom, double e, double n)
{
    zoom_ = std::max(1.0, std::min(kMaxZoom, zoom));
    /* keep the instrumented range in view: the view centre stays within it */
    const double r = std::hypot(e, n), rMax = maxRange_m_ * (1.0 - 1.0 / zoom_);
    if (r > rMax) {
        e *= rMax / r;
        n *= rMax / r;
    }
    viewE_ = e;
    viewN_ = n;
    backgroundDirty_ = true;
    update();
}

void PpiWidget::zoomAt(double factor, QPointF screen)
{
    layoutDisplay();
    const QPointF g = toGround(screen);          /* this ground point stays under the cursor */
    const double z = std::max(1.0, std::min(kMaxZoom, zoom_ * factor));
    const double k = radius_px_ * z / maxRange_m_;
    setView(z, g.x() - (screen.x() - centre_.x()) / k, g.y() + (screen.y() - centre_.y()) / k);
}

void PpiWidget::panPixels(QPointF delta)
{
    layoutDisplay();
    const double k = scale();
    setView(zoom_, viewE_ - delta.x() / k, viewN_ + delta.y() / k);
}

void PpiWidget::resetView()
{
    setView(1.0, 0.0, 0.0);
}

/* ---- events -------------------------------------------------------------- */

void PpiWidget::closeEvent(QCloseEvent *event)
{
    if (!finished_ && stopFlag_) *stopFlag_ = 1;   /* ask the worker to stop */
    event->accept();
}

void PpiWidget::wheelEvent(QWheelEvent *event)
{
    const double steps = event->angleDelta().y() / 120.0;
    if (steps != 0.0) zoomAt(std::pow(1.25, steps), event_pos(event));
    event->accept();
}

void PpiWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() != Qt::LeftButton) return;
    dragging_ = true;
    dragMoved_ = false;
    dragStart_ = event_pos(event);
    dragE_ = viewE_;
    dragN_ = viewN_;
}

void PpiWidget::mouseMoveEvent(QMouseEvent *event)
{
    const QPointF pos = event_pos(event);
    if (dragging_) {
        const QPointF d = pos - dragStart_;
        if (std::fabs(d.x()) > 3.0 || std::fabs(d.y()) > 3.0) dragMoved_ = true;
        if (dragMoved_) {
            const double k = scale();
            setView(zoom_, dragE_ - d.x() / k, dragN_ + d.y() / k);
            setCursor(Qt::ClosedHandCursor);
        }
    }
    updateCursor(pos);
    const int i = dragging_ ? -1 : pick(pos);
    const uint64_t id = i < 0 ? 0 : plots_[(size_t)i].id;
    if (id != hoverId_) {
        hoverId_ = id;
        if (id) QToolTip::showText(mapToGlobal(pos.toPoint()), plot_text(plots_[(size_t)i].det), this);
        else QToolTip::hideText();
        update();
    }
}

void PpiWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() != Qt::LeftButton || !dragging_) return;
    if (!dragMoved_) clickAt(event_pos(event));
    dragging_ = false;
    unsetCursor();
}

void PpiWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) resetView();
}

void PpiWidget::leaveEvent(QEvent *event)
{
    (void)event;
    cursorText_.clear();
    hoverId_ = 0;
    update();
}

void PpiWidget::keyPressEvent(QKeyEvent *event)
{
    const QPointF mid = centre_;
    const double step = 0.1 * std::min(width(), height());
    switch (event->key()) {
    case Qt::Key_Plus: case Qt::Key_Equal: zoomAt(1.5, mid); break;
    case Qt::Key_Minus:                    zoomAt(1.0 / 1.5, mid); break;
    case Qt::Key_Left:  panPixels(QPointF(step, 0.0)); break;
    case Qt::Key_Right: panPixels(QPointF(-step, 0.0)); break;
    case Qt::Key_Up:    panPixels(QPointF(0.0, step)); break;
    case Qt::Key_Down:  panPixels(QPointF(0.0, -step)); break;
    case Qt::Key_R: case Qt::Key_Home: resetView(); break;
    case Qt::Key_C:     pinned_.clear(); update(); break;
    case Qt::Key_Space: setPaused(!paused_); break;
    case Qt::Key_S: {
        const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("Save the PPI"),
                                                          QStringLiteral("ppi.png"), QStringLiteral("PNG (*.png)"));
        if (!path.isEmpty()) saveSnapshot(path);
        break;
    }
    default: QWidget::keyPressEvent(event);
    }
}

void PpiWidget::contextMenuEvent(QContextMenuEvent *event)
{
    QMenu menu(this);
    auto toggle = [&menu](const QString &text, bool on, std::function<void(bool)> f) {
        QAction *a = menu.addAction(text);
        a->setCheckable(true);
        a->setChecked(on);
        QObject::connect(a, &QAction::toggled, f);
    };
    toggle(QStringLiteral("Short-pulse plots"), showShort_, [this](bool on) { setShowShort(on); });
    toggle(QStringLiteral("Long-pulse plots"), showLong_, [this](bool on) { setShowLong(on); });
    toggle(QStringLiteral("Plot labels (range, velocity)"), showLabels_, [this](bool on) { showLabels_ = on; update(); });
    toggle(QStringLiteral("Range rings and spokes"), showRings_,
           [this](bool on) { showRings_ = on; backgroundDirty_ = true; update(); });
    toggle(QStringLiteral("Antenna trail"), showTrail_, [this](bool on) { showTrail_ = on; update(); });

    QMenu *hist = menu.addMenu(QStringLiteral("Plots shown"));
    QActionGroup *histGroup = new QActionGroup(hist);
    const std::pair<const char *, History> histItems[] = {
        { "Until the antenna passes again", History::UntilSwept },
        { "Two scans", History::TwoScans },
        { "All (until cleared)", History::All } };
    for (const auto &item : histItems) {
        QAction *a = hist->addAction(QString::fromLatin1(item.first));
        a->setCheckable(true);
        a->setChecked(history_ == item.second);
        histGroup->addAction(a);
        const History h = item.second;
        QObject::connect(a, &QAction::triggered, [this, h]() { setHistory(h); });
    }
    QMenu *snr = menu.addMenu(QStringLiteral("Minimum SNR"));
    QActionGroup *snrGroup = new QActionGroup(snr);
    for (double dB : { -1e9, 10.0, 13.0, 15.0, 20.0, 25.0, 30.0 }) {
        QAction *a = snr->addAction(dB < -1e8 ? QStringLiteral("All") : QStringLiteral("%1 dB").arg(dB));
        a->setCheckable(true);
        a->setChecked(minSnr_dB_ == dB);
        snrGroup->addAction(a);
        QObject::connect(a, &QAction::triggered, [this, dB]() { setMinSnr(dB); });
    }
    menu.addSeparator();
    QObject::connect(menu.addAction(QStringLiteral("Reset view (double click)")), &QAction::triggered,
                     [this]() { resetView(); });
    QObject::connect(menu.addAction(QStringLiteral("Clear plots")), &QAction::triggered, [this]() { clearPlots(); });
    QObject::connect(menu.addAction(QStringLiteral("Remove data tips (C)")), &QAction::triggered,
                     [this]() { pinned_.clear(); update(); });
    QObject::connect(menu.addAction(paused_ ? QStringLiteral("Resume display (Space)") : QStringLiteral("Pause display (Space)")),
                     &QAction::triggered, [this]() { setPaused(!paused_); });
    QObject::connect(menu.addAction(QStringLiteral("Save as PNG ... (S)")), &QAction::triggered, [this]() {
        const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("Save the PPI"),
                                                          QStringLiteral("ppi.png"), QStringLiteral("PNG (*.png)"));
        if (!path.isEmpty()) saveSnapshot(path);
    });
    menu.exec(event->globalPos());
}

void PpiWidget::updateCursor(QPointF screen)
{
    if (!viewRect_.contains(screen)) {
        cursorText_.clear();
    } else {
        const QPointF g = toGround(screen);
        cursorText_ = QStringLiteral("Cursor %1, %2   (E %3, N %4)")
                          .arg(km_text(std::hypot(g.x(), g.y())))
                          .arg(deg_text(wrap360(std::atan2(g.x(), g.y()) / kDeg), 1))
                          .arg(km_text(g.x()))
                          .arg(km_text(g.y()));
    }
    update();
}

/* ---- drawing ------------------------------------------------------------- */

/* Square polar display below the status text, with room for the labels. */
void PpiWidget::layoutDisplay()
{
    const double w = width(), h = height();
    const double top = kHeaderPx + kLabelPx;
    const double side = std::max(50.0, std::min(w - 2.0 * kLabelPx, h - top - kLabelPx));
    radius_px_ = side / 2.0;
    centre_ = QPointF(w / 2.0, top + (h - top - kLabelPx) / 2.0);
    viewRect_ = QRectF(0.0, kHeaderPx, w, h - kHeaderPx);
}

/* Farthest point of the ray from the radar along az_deg (up to rEnd_px) that
 * lies inside rect; false if the ray misses the rectangle. */
bool ray_end_in_rect(QPointF origin, double az_deg, double rEnd_px, const QRectF &rect, QPointF *end)
{
    const double dx = std::sin(az_deg * kDeg), dy = -std::cos(az_deg * kDeg);
    double t0 = 0.0, t1 = rEnd_px;                         /* Liang-Barsky clipping */
    const double p[4] = { -dx, dx, -dy, dy };
    const double q[4] = { origin.x() - rect.left(), rect.right() - origin.x(),
                          origin.y() - rect.top(), rect.bottom() - origin.y() };
    for (int i = 0; i < 4; ++i) {
        if (std::fabs(p[i]) < 1e-12) {
            if (q[i] < 0.0) return false;
            continue;
        }
        const double t = q[i] / p[i];
        if (p[i] < 0.0) t0 = std::max(t0, t);
        else t1 = std::min(t1, t);
    }
    if (t0 > t1) return false;
    *end = QPointF(origin.x() + t1 * dx, origin.y() + t1 * dy);
    return true;
}

void PpiWidget::drawBackground()
{
    const qreal dpr = devicePixelRatioF();
    background_ = QPixmap(size() * dpr);
    background_.setDevicePixelRatio(dpr);
    background_.fill(Qt::black);
    backgroundSize_ = size();
    backgroundDirty_ = false;

    QPainter p(&background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setClipRect(viewRect_);
    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    p.setFont(font);
    const double k = scale(), rim_px = maxRange_m_ * k;
    const QPointF origin = toScreen(0.0, 0.0);

    if (showRings_) {
        /* azimuth spokes every 30 deg; labels outside the rim or, when zoomed,
         * at the edge of the view */
        const QRectF inner = viewRect_.adjusted(14.0, 10.0, -14.0, -10.0);
        for (int az = 0; az < 360; az += 30) {
            QPointF end;
            if (!ray_end_in_rect(origin, az, rim_px, viewRect_, &end)) continue;
            p.setPen(QPen(kColorGrid, 1.0));
            p.drawLine(origin, polar(az, maxRange_m_));
            QPointF c;
            if (zoom_ <= 1.0 + 1e-9) c = polar(az, maxRange_m_ * (1.0 + 15.0 / rim_px));
            else if (!ray_end_in_rect(origin, az, rim_px, inner, &c)) continue;
            const QString text = deg_text(az, 0);
            p.setPen(kColorLabel);
            p.drawText(QPointF(c.x() - fm.horizontalAdvance(text) / 2.0, c.y() + fm.ascent() / 2.0 - 1.0), text);
        }
        /* range rings for the visible range, labels between the 0 and 30 deg spokes */
        const double step = ring_step_m(maxRange_m_ / zoom_);
        for (double r = step; r <= maxRange_m_ * (1.0 + 1e-9); r += step) {
            p.setPen(QPen(kColorGrid, 1.0));
            p.setBrush(Qt::NoBrush);
            p.drawEllipse(origin, r * k, r * k);
            p.setPen(kColorLabel);
            p.drawText(polar(15.0, r) + QPointF(3.0, -3.0), km_text(r));
        }
    }
    /* rim of the instrumented range */
    p.setPen(QPen(kColorRim, 1.5));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(origin, rim_px, rim_px);
}

/* Antenna line with a fading trail behind it (the antenna turns clockwise). */
void PpiWidget::drawAntenna(QPainter &p) const
{
    if (!antenna_.valid) return;
    const QPointF origin = toScreen(0.0, 0.0);
    const double rim_px = maxRange_m_ * scale();
    if (showTrail_) {
        const QRectF disc(origin.x() - rim_px, origin.y() - rim_px, 2.0 * rim_px, 2.0 * rim_px);
        const double qtAngle = 90.0 - antenna_.az_deg;   /* Qt: counter-clockwise from East */
        QConicalGradient trail(origin, qtAngle);
        QColor head(kColorAntenna), tail(kColorAntenna);
        head.setAlpha(90);
        tail.setAlpha(0);
        trail.setColorAt(0.0, head);
        trail.setColorAt(kTrailDeg / 360.0, tail);
        trail.setColorAt(1.0, tail);
        p.setPen(Qt::NoPen);
        p.setBrush(trail);
        p.drawPie(disc, (int)std::lround(qtAngle * 16.0), (int)std::lround(kTrailDeg * 16.0));
    }
    p.setPen(QPen(kColorAntenna, 2.0));
    p.drawLine(origin, polar(antenna_.az_deg, maxRange_m_));
}

/* Filled circles; older plots are dimmer (phosphor-like decay over a scan). */
void PpiWidget::drawPlots(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 11);
    p.setFont(font);
    for (const PlotMark &m : plots_) {
        if (!shown(m)) continue;
        const QPointF s = toScreen(m.east_m, m.north_m);
        if (!viewRect_.adjusted(-10, -10, 10, 10).contains(s)) continue;
        QColor c = (m.pulse == PULSE_LONG) ? kColorLong : kColorShort;
        const double age = std::min(1.0, turned(m) / 360.0);
        c.setAlpha(255 - (int)(130.0 * age));
        p.setPen(Qt::NoPen);
        p.setBrush(c);
        p.drawEllipse(s, kPlotRadius, kPlotRadius);
        if (m.id == hoverId_) {                         /* highlight the plot under the cursor */
            p.setPen(QPen(Qt::white, 1.5));
            p.setBrush(Qt::NoBrush);
            p.drawEllipse(s, kPlotRadius + 4.0, kPlotRadius + 4.0);
        }
        if (showLabels_) {
            p.setPen(kColorLabel);
            p.drawText(s + QPointF(7.0, 4.0), QStringLiteral("%1 km, %2 m/s")
                                                  .arg(m.det.range_m / 1e3, 0, 'f', 1)
                                                  .arg(m.det.vr_mps, 0, 'f', 1));
        }
    }
}

/* Pinned data tips: a box with the plot values next to the plot. */
void PpiWidget::drawPins(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 11);
    const QFontMetrics fm(font);
    p.setFont(font);
    for (const PlotMark &m : plots_) {
        if (!pinned_.count(m.id)) continue;
        const QPointF s = toScreen(m.east_m, m.north_m);
        QString text = plot_text(m.det);
        text.replace(QStringLiteral("<b>"), QString()).replace(QStringLiteral("</b>"), QString());
        const QStringList rows = text.split(QStringLiteral("<br>"));
        int w = 0;
        for (const QString &r : rows) w = std::max(w, fm.horizontalAdvance(r));
        const double rowH = fm.height(), boxW = w + 12.0, boxH = rows.size() * rowH + 8.0;
        QRectF box(s.x() + 12.0, s.y() + 8.0, boxW, boxH);
        if (box.right() > width()) box.moveRight(s.x() - 12.0);          /* stay inside the widget */
        if (box.bottom() > height()) box.moveBottom(s.y() - 8.0);
        p.setPen(QPen(kColorPin, 1.5));
        p.setBrush(Qt::NoBrush);
        p.drawRect(QRectF(s.x() - 7.0, s.y() - 7.0, 14.0, 14.0));
        p.drawLine(s, box.topLeft());
        p.setBrush(QColor(10, 40, 10, 225));
        p.drawRect(box);
        p.setPen(Qt::white);
        for (int i = 0; i < rows.size(); ++i)
            p.drawText(QPointF(box.left() + 6.0, box.top() + 4.0 + fm.ascent() + i * rowH), rows[i]);
    }
}

void PpiWidget::drawStatus(QPainter &p) const
{
    QString line1;
    if (antenna_.valid)
        line1 = QStringLiteral("Scan %1    t = %2 s    Antenna az %3    Plots %4 / %5")
                    .arg(antenna_.scan)
                    .arg(antenna_.time_s, 0, 'f', 2)
                    .arg(deg_text(antenna_.az_deg, 1))
                    .arg((qulonglong)visibleCount())
                    .arg((qulonglong)plots_.size());
    else
        line1 = finished_ ? QStringLiteral("No antenna data received")
                          : QStringLiteral("Waiting for data ...");
    if (paused_) line1 += QStringLiteral("    [PAUSED]");

    QString line2 = QStringLiteral("Range %1, zoom x%2, rings every %3    ")
                        .arg(km_text(maxRange_m_))
                        .arg(zoom_, 0, 'f', 1)
                        .arg(km_text(ring_step_m(maxRange_m_ / zoom_)));
    if (!cursorText_.isEmpty())
        line2 += cursorText_;
    else if (!finished_)
        line2 += QStringLiteral("receiving");
    else if (workerCode_ == 0)
        line2 += QStringLiteral("processing finished");
    else
        line2 += QStringLiteral("processing stopped with error %1").arg(workerCode_);

    p.setPen(Qt::white);
    p.setFont(pixel_font(font(), 16, true));
    p.drawText(QPointF(12.0, 22.0), line1);
    p.setPen(kColorLabel);
    p.setFont(pixel_font(font(), 12));
    p.drawText(QPointF(12.0, 41.0), line2);
}

/* Small legend in the bottom-left corner, with the mouse help. */
void PpiWidget::drawLegend(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    const QString names[4] = { QStringLiteral("Plot, short pulse"), QStringLiteral("Plot, long pulse"),
                               QStringLiteral("Antenna"),
                               QStringLiteral("Wheel zoom, drag pan, right click options") };
    int textWidth = 0;
    for (const QString &n : names) textWidth = std::max(textWidth, fm.horizontalAdvance(n));

    const double rowH = 18.0, pad = 8.0;
    const QRectF box(8.0, height() - 8.0 - 4.0 * rowH - 2.0 * pad, 2.0 * pad + 26.0 + textWidth, 4.0 * rowH + 2.0 * pad);
    p.setPen(QPen(kColorGrid, 1.0));
    p.setBrush(QColor(20, 20, 20, 220));
    p.drawRect(box);

    p.setFont(font);
    for (int i = 0; i < 4; ++i) {
        const double y = box.top() + pad + rowH * (i + 0.5);
        const QPointF mark(box.left() + pad + 8.0, y);
        if (i < 2) {
            p.setPen(Qt::NoPen);
            p.setBrush(i == 0 ? kColorShort : kColorLong);
            p.drawEllipse(mark, kPlotRadius, kPlotRadius);
        } else if (i == 2) {
            p.setPen(QPen(kColorAntenna, 2.0));
            p.drawLine(mark - QPointF(8.0, 0.0), mark + QPointF(8.0, 0.0));
        }
        p.setPen(i == 3 ? kColorLabel : QColor(Qt::white));
        p.drawText(QPointF(box.left() + pad + (i == 3 ? 0.0 : 26.0), y + fm.ascent() / 2.0 - 1.0), names[i]);
    }
}

void PpiWidget::paintEvent(QPaintEvent *event)
{
    (void)event;
    layoutDisplay();
    if (background_.isNull() || backgroundDirty_ || backgroundSize_ != size()) drawBackground();

    QPainter p(this);
    p.drawPixmap(0, 0, background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.save();
    p.setClipRect(viewRect_);
    drawAntenna(p);
    drawPlots(p);
    p.restore();
    drawPins(p);
    drawStatus(p);
    drawLegend(p);
}

bool save_snapshot(PpiWidget &widget, const char *path)
{
    return widget.saveSnapshot(QString::fromLocal8Bit(path));
}

/* Without a display (e.g. a Linux server) a snapshot is rendered off screen. */
void use_offscreen_without_display(const char *snapshotPath)
{
#if defined(__linux__)
    if (snapshotPath && qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM") && qEnvironmentVariableIsEmpty("DISPLAY")
        && qEnvironmentVariableIsEmpty("WAYLAND_DISPLAY"))
        qputenv("QT_QPA_PLATFORM", "offscreen");
#else
    (void)snapshotPath;
#endif
}

} /* namespace */

/* ------------------------------------------------------------------------ */
/* Public C interface                                                        */
/* ------------------------------------------------------------------------ */

extern "C" void ppi_post_antenna(double az_deg, double time_s, uint32_t scan)
{
    if (!std::isfinite(az_deg)) return;
    az_deg = wrap360(az_deg);
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    AntennaState &a = s.antenna;
    if (a.valid && time_s < a.time_s - 1.0) {      /* time went back: a new run */
        a = AntennaState();
        s.newPlots.clear();
        s.clear = true;
    }
    if (a.valid) a.sweep_deg += antenna_step(a, az_deg, scan);
    a.valid = true;
    a.az_deg = az_deg;
    a.time_s = time_s;
    a.scan = scan;
}

extern "C" void ppi_post_plot(const Detection *plot)
{
    if (!plot || !std::isfinite(plot->range_m) || !std::isfinite(plot->az_deg)) return;
    const double az = wrap360(plot->az_deg);
    double ground_m = plot->range_m;
    if (std::isfinite(plot->el_deg)) ground_m *= std::cos(plot->el_deg * kDeg);

    PlotMark m;
    m.east_m = ground_m * std::sin(az * kDeg);
    m.north_m = ground_m * std::cos(az * kDeg);
    m.pulse = (plot->pulse == PULSE_LONG) ? PULSE_LONG : PULSE_SHORT;
    m.det = *plot;

    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    const AntennaState &a = s.antenna;
    /* the plot azimuth is close to (just behind) the antenna */
    m.sweep_deg = a.sweep_deg + (a.valid ? wrap180(az - a.az_deg) : 0.0);
    m.scan = a.scan;
    m.id = s.nextId++;
    if (s.newPlots.size() < kMaxPending) s.newPlots.push_back(m);
}

extern "C" int ppi_run(int argc, char **argv, double maxRange_m, int (*worker)(void *), void *arg,
                       volatile sig_atomic_t *stopFlag, const char *snapshotPath)
{
    static char defaultName[] = "radar_receiver";
    static char *defaultArgv[] = { defaultName, nullptr };
    int qtArgc = argc;                  /* must outlive the QApplication */
    char **qtArgv = argv;
    if (qtArgc < 1 || !qtArgv || !qtArgv[0]) {
        qtArgc = 1;
        qtArgv = defaultArgv;
    }

    reset_shared();
    use_offscreen_without_display(snapshotPath);
    QApplication app(qtArgc, qtArgv);

    PpiWidget widget(maxRange_m, stopFlag);
    widget.setWindowTitle(QStringLiteral("Radar PPI"));
    widget.resize(800, 800);
    widget.show();

    std::atomic<bool> workerDone(false);
    int workerCode = 0;
    std::thread thread([&]() {
        if (worker) workerCode = worker(arg);
        workerDone.store(true);
    });

    bool finishSeen = false;
    QTimer timer;
    QObject::connect(&timer, &QTimer::timeout, [&]() {
        widget.refresh();
        if (finishSeen || !workerDone.load()) return;
        finishSeen = true;
        if (snapshotPath) {             /* the snapshot shows all data */
            widget.setPaused(false);
            widget.refresh();
        }
        widget.setFinished(workerCode);
        if (snapshotPath)
            QCoreApplication::exit(save_snapshot(widget, snapshotPath) ? 0 : 1);
        else if (stopFlag && *stopFlag)
            QCoreApplication::quit();   /* stopped from outside (e.g. Ctrl+C) */
        /* otherwise the window stays open until the user closes it */
    });
    timer.start(kRefreshMs);

    int code = app.exec();
    timer.stop();
    if (!workerDone.load() && stopFlag) *stopFlag = 1;   /* window closed first */
    thread.join();
    if (code == 0) code = workerCode;   /* keep the error code of the worker */
    return code;
}
