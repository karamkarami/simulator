/*
 * ppi_qt.cpp - Qt PPI (plan position indicator) of the radar receiver.
 *
 * The worker thread posts antenna positions and plots with ppi_post_antenna()
 * and ppi_post_plot(). These calls only store the data under a mutex; a timer
 * in the GUI thread (25 Hz) takes the new data and repaints the widget.
 *
 * Plot persistence: a plot stays on the screen until the antenna sweeps over
 * its azimuth again in a later scan, like on a real PPI. For this the antenna
 * rotation is accumulated ("sweep", degrees turned since the first update), so
 * the scan boundary (the start azimuth of the antenna) does not matter.
 *
 * The file does not use Q_OBJECT, so it builds with or without moc, with
 * Qt 5 (>= 5.12) and Qt 6.
 */
#include "radar_rx.h"

#include <QApplication>
#include <QCloseEvent>
#include <QColor>
#include <QConicalGradient>
#include <QFont>
#include <QFontMetrics>
#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QPixmap>
#include <QPointF>
#include <QRectF>
#include <QString>
#include <QTimer>
#include <QWidget>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <deque>
#include <mutex>
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
constexpr double kTrailDeg    = 20.0;     /* fading trail behind the antenna line   */

const QColor kColorShort(0, 230, 255);    /* cyan                                   */
const QColor kColorLong(255, 230, 0);     /* yellow                                 */
const QColor kColorAntenna(0, 220, 0);    /* green                                  */
const QColor kColorGrid(70, 70, 70);
const QColor kColorRim(120, 120, 120);
const QColor kColorLabel(150, 150, 150);

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
    double   east_m = 0.0;      /* ground position                              */
    double   north_m = 0.0;
    int      pulse = PULSE_SHORT;
    double   sweep_deg = 0.0;   /* antenna sweep at which it pointed at the plot */
    uint32_t scan = 0;          /* antenna scan when the plot arrived           */
};

struct Shared {
    std::mutex            mutex;
    AntennaState          antenna;
    std::vector<PlotMark> newPlots;   /* plots not yet taken by the GUI thread */
    bool                  clear = false;   /* a new run started: clear the display */
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

/* Ring spacing: 10 km, 20 km above 60 km (larger/smaller for extreme ranges). */
double ring_step_m(double maxRange_m)
{
    const double km = maxRange_m / 1e3;
    if (km > 500.0) return 100e3;
    if (km > 200.0) return 50e3;
    if (km > 60.0)  return 20e3;
    if (km > 30.0)  return 10e3;
    if (km > 15.0)  return 5e3;
    if (km > 6.0)   return 2e3;
    return 1e3;
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

/* ------------------------------------------------------------------------ */
/* PPI widget                                                                */
/* ------------------------------------------------------------------------ */

class PpiWidget : public QWidget {
public:
    PpiWidget(double maxRange_m, volatile sig_atomic_t *stopFlag)
        : maxRange_m_(maxRange_m), ringStep_m_(ring_step_m(maxRange_m)), stopFlag_(stopFlag)
    {
        setAttribute(Qt::WA_OpaquePaintEvent);
        setMinimumSize(300, 300);
    }

    void refresh();                     /* takes the posted data (GUI thread) */
    void setFinished(int workerCode)
    {
        finished_ = true;
        workerCode_ = workerCode;
        update();
    }

protected:
    void paintEvent(QPaintEvent *event) override;
    void closeEvent(QCloseEvent *event) override;

private:
    void   layoutDisplay();
    void   drawBackground();
    void   drawAntenna(QPainter &p) const;
    void   drawPlots(QPainter &p) const;
    void   drawStatus(QPainter &p) const;
    void   drawLegend(QPainter &p) const;
    bool   expired(const PlotMark &m) const;
    double turned(const PlotMark &m) const { return std::fabs(antenna_.sweep_deg - m.sweep_deg); }
    QPointF toScreen(double east_m, double north_m) const
    {
        const double k = radius_px_ / maxRange_m_;
        return QPointF(centre_.x() + k * east_m, centre_.y() - k * north_m);
    }
    QPointF polar(double az_deg, double r_px) const
    {
        return QPointF(centre_.x() + r_px * std::sin(az_deg * kDeg), centre_.y() - r_px * std::cos(az_deg * kDeg));
    }

    double                maxRange_m_;
    double                ringStep_m_;
    volatile sig_atomic_t *stopFlag_;
    AntennaState          antenna_;
    std::deque<PlotMark>  plots_;
    std::vector<PlotMark> incoming_;    /* buffer swapped with Shared::newPlots */
    bool                  finished_ = false;
    int                   workerCode_ = 0;
    QPointF               centre_;
    double                radius_px_ = 1.0;
    QPixmap               background_;  /* rings, spokes and labels            */
    QSize                 backgroundSize_;
};

void PpiWidget::refresh()
{
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
    if (clear) plots_.clear();
    plots_.insert(plots_.end(), incoming_.begin(), incoming_.end());
    plots_.erase(std::remove_if(plots_.begin(), plots_.end(),
                                [this](const PlotMark &m) { return expired(m); }),
                 plots_.end());
    while (plots_.size() > kMaxPlots) plots_.pop_front();
    update();
}

/* A plot disappears when the antenna has turned a full revolution since it
 * pointed at the plot, i.e. it sweeps over the plot azimuth in the next scan.
 * Safety net for gaps in the data: two scans later it is gone in any case. */
bool PpiWidget::expired(const PlotMark &m) const
{
    const long long scans = (long long)antenna_.scan - (long long)m.scan;
    return turned(m) >= 360.0 || scans >= 2 || scans <= -2;
}

void PpiWidget::closeEvent(QCloseEvent *event)
{
    if (!finished_ && stopFlag_) *stopFlag_ = 1;   /* ask the worker to stop */
    event->accept();
}

/* Square polar display below the status text, with room for the labels. */
void PpiWidget::layoutDisplay()
{
    const double w = width(), h = height();
    const double top = kHeaderPx + kLabelPx;
    const double side = std::max(50.0, std::min(w - 2.0 * kLabelPx, h - top - kLabelPx));
    radius_px_ = side / 2.0;
    centre_ = QPointF(w / 2.0, top + (h - top - kLabelPx) / 2.0);
    if (!(maxRange_m_ > 0.0)) maxRange_m_ = 60e3;   /* also catches NaN */
}

void PpiWidget::drawBackground()
{
    const qreal dpr = devicePixelRatioF();
    background_ = QPixmap(size() * dpr);
    background_.setDevicePixelRatio(dpr);
    background_.fill(Qt::black);
    backgroundSize_ = size();

    QPainter p(&background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    p.setFont(font);

    /* azimuth spokes every 30 deg with labels outside the rim */
    for (int az = 0; az < 360; az += 30) {
        p.setPen(QPen(kColorGrid, 1.0));
        p.drawLine(centre_, polar(az, radius_px_));
        const QString text = deg_text(az, 0);
        const QPointF c = polar(az, radius_px_ + 15.0);
        p.setPen(kColorLabel);
        p.drawText(QPointF(c.x() - fm.horizontalAdvance(text) / 2.0, c.y() + fm.ascent() / 2.0 - 1.0), text);
    }

    /* range rings with labels between the 0 and 30 deg spokes */
    const double k = radius_px_ / maxRange_m_;
    for (double r = ringStep_m_; r <= maxRange_m_ * (1.0 + 1e-9); r += ringStep_m_) {
        p.setPen(QPen(kColorGrid, 1.0));
        p.setBrush(Qt::NoBrush);
        p.drawEllipse(centre_, r * k, r * k);
        p.setPen(kColorLabel);
        p.drawText(polar(15.0, r * k) + QPointF(3.0, -3.0), km_text(r));
    }

    /* rim of the instrumented range */
    p.setPen(QPen(kColorRim, 1.5));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(centre_, radius_px_, radius_px_);
}

/* Antenna line with a fading trail behind it (the antenna turns clockwise). */
void PpiWidget::drawAntenna(QPainter &p) const
{
    if (!antenna_.valid) return;
    const QRectF disc(centre_.x() - radius_px_, centre_.y() - radius_px_, 2.0 * radius_px_, 2.0 * radius_px_);
    const double qtAngle = 90.0 - antenna_.az_deg;   /* Qt: counter-clockwise from East */
    QConicalGradient trail(centre_, qtAngle);
    QColor head(kColorAntenna), tail(kColorAntenna);
    head.setAlpha(90);
    tail.setAlpha(0);
    trail.setColorAt(0.0, head);
    trail.setColorAt(kTrailDeg / 360.0, tail);
    trail.setColorAt(1.0, tail);
    p.setPen(Qt::NoPen);
    p.setBrush(trail);
    p.drawPie(disc, (int)std::lround(qtAngle * 16.0), (int)std::lround(kTrailDeg * 16.0));

    p.setPen(QPen(kColorAntenna, 2.0));
    p.drawLine(centre_, polar(antenna_.az_deg, radius_px_));
}

/* Filled circles; older plots are dimmer (phosphor-like decay over a scan). */
void PpiWidget::drawPlots(QPainter &p) const
{
    p.setPen(Qt::NoPen);
    for (const PlotMark &m : plots_) {
        QColor c = (m.pulse == PULSE_LONG) ? kColorLong : kColorShort;
        const double age = std::min(1.0, turned(m) / 360.0);
        c.setAlpha(255 - (int)(130.0 * age));
        p.setBrush(c);
        p.drawEllipse(toScreen(m.east_m, m.north_m), kPlotRadius, kPlotRadius);
    }
}

void PpiWidget::drawStatus(QPainter &p) const
{
    QString line1;
    if (antenna_.valid)
        line1 = QStringLiteral("Scan %1    t = %2 s    Antenna az %3    Plots %4")
                    .arg(antenna_.scan)
                    .arg(antenna_.time_s, 0, 'f', 2)
                    .arg(deg_text(antenna_.az_deg, 1))
                    .arg((qulonglong)plots_.size());
    else
        line1 = finished_ ? QStringLiteral("No antenna data received")
                          : QStringLiteral("Waiting for data ...");

    QString line2 = QStringLiteral("Range 0 - %1, rings every %2    ").arg(km_text(maxRange_m_), km_text(ringStep_m_));
    if (!finished_)
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

/* Small legend in the bottom-left corner. */
void PpiWidget::drawLegend(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    const QString names[3] = { QStringLiteral("Plot, short pulse"), QStringLiteral("Plot, long pulse"),
                               QStringLiteral("Antenna") };
    int textWidth = 0;
    for (const QString &n : names) textWidth = std::max(textWidth, fm.horizontalAdvance(n));

    const double rowH = 18.0, pad = 8.0;
    const QRectF box(8.0, height() - 8.0 - 3.0 * rowH - 2.0 * pad, 2.0 * pad + 26.0 + textWidth, 3.0 * rowH + 2.0 * pad);
    p.setPen(QPen(kColorGrid, 1.0));
    p.setBrush(QColor(20, 20, 20, 220));
    p.drawRect(box);

    p.setFont(font);
    for (int i = 0; i < 3; ++i) {
        const double y = box.top() + pad + rowH * (i + 0.5);
        const QPointF mark(box.left() + pad + 8.0, y);
        if (i < 2) {
            p.setPen(Qt::NoPen);
            p.setBrush(i == 0 ? kColorShort : kColorLong);
            p.drawEllipse(mark, kPlotRadius, kPlotRadius);
        } else {
            p.setPen(QPen(kColorAntenna, 2.0));
            p.drawLine(mark - QPointF(8.0, 0.0), mark + QPointF(8.0, 0.0));
        }
        p.setPen(Qt::white);
        p.drawText(QPointF(box.left() + pad + 26.0, y + fm.ascent() / 2.0 - 1.0), names[i]);
    }
}

void PpiWidget::paintEvent(QPaintEvent *event)
{
    (void)event;
    layoutDisplay();
    if (background_.isNull() || backgroundSize_ != size()) drawBackground();

    QPainter p(this);
    p.drawPixmap(0, 0, background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    drawAntenna(p);
    drawPlots(p);
    drawStatus(p);
    drawLegend(p);
}

bool save_snapshot(QWidget &widget, const char *path)
{
    const QPixmap pixmap = widget.grab();
    if (pixmap.save(QString::fromLocal8Bit(path), "PNG")) return true;
    std::fprintf(stderr, "PPI: cannot save the snapshot to %s\n", path);
    return false;
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

    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    const AntennaState &a = s.antenna;
    /* the plot azimuth is close to (just behind) the antenna */
    m.sweep_deg = a.sweep_deg + (a.valid ? wrap180(az - a.az_deg) : 0.0);
    m.scan = a.scan;
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
