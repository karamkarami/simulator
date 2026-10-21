/*
 * test_ppi_interactive.cpp - interaction with the Qt PPI.
 *
 * Includes ppi_qt.cpp to reach the widget, sends real mouse, wheel and key
 * events (off-screen) and checks: hover data tips, cursor readout, wheel zoom
 * around the cursor, drag pan, click to pin a data tip, double click reset,
 * keys, show/hide and SNR filters, plot persistence options, pause. Saves
 * snapshots (argv[1]: folder, optional) of the zoomed PPI with a pinned tip.
 */
#include "../src/ppi_qt.cpp"

#include <cstring>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char *what)
{
    std::printf("  %-60s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) failures++;
}

void send_mouse(QWidget &w, QEvent::Type type, QPointF pos, Qt::MouseButton button)
{
    const Qt::MouseButtons buttons = (type == QEvent::MouseButtonRelease || type == QEvent::MouseMove)
                                         ? (type == QEvent::MouseMove && button != Qt::NoButton ? Qt::MouseButtons(button) : Qt::MouseButtons(Qt::NoButton))
                                         : Qt::MouseButtons(button);
    QMouseEvent e(type, pos, QPointF(w.mapToGlobal(pos.toPoint())), button, buttons, Qt::NoModifier);
    QCoreApplication::sendEvent(&w, &e);
}

void send_wheel(QWidget &w, QPointF pos, int angle)
{
    QWheelEvent e(pos, w.mapToGlobal(pos.toPoint()), QPoint(), QPoint(0, angle), Qt::NoButton, Qt::NoModifier,
                  Qt::NoScrollPhase, false);
    QCoreApplication::sendEvent(&w, &e);
}

void send_key(QWidget &w, int key, const char *text)
{
    QKeyEvent e(QEvent::KeyPress, key, Qt::NoModifier, QString::fromLatin1(text));
    QCoreApplication::sendEvent(&w, &e);
}

Detection make_plot(double range_m, double az_deg, double el_deg, double vr, double snr, int pulse, uint32_t number)
{
    Detection d;
    std::memset(&d, 0, sizeof d);
    d.range_m = range_m;  d.az_deg = az_deg;  d.el_deg = el_deg;  d.vr_mps = vr;  d.snr_dB = snr;
    d.pulse = pulse;  d.hits = 70;  d.beam = 4;  d.rIdx = 123;  d.dIdx = 45;  d.time_s = 1.5;
    d.pulseNumber = number;
    return d;
}

QPointF ground(const Detection &d)
{
    const double g = d.range_m * (std::isfinite(d.el_deg) ? std::cos(d.el_deg * kDeg) : 1.0);
    return QPointF(g * std::sin(d.az_deg * kDeg), g * std::cos(d.az_deg * kDeg));
}

} /* namespace */

int main(int argc, char **argv)
{
    if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM")) qputenv("QT_QPA_PLATFORM", "offscreen");
    const std::string folder = argc > 1 ? argv[1] : ".";
    QApplication app(argc, argv);
    reset_shared();
    PpiWidget w(60e3, nullptr);
    w.resize(800, 800);
    w.show();

    const Detection a = make_plot(10e3, 30.0, NAN, 12.5, 25.0, PULSE_SHORT, 101);
    const Detection b = make_plot(40e3, 200.0, 3.0, -7.25, 14.0, PULSE_LONG, 102);
    const Detection c = make_plot(30e3, 95.0, 5.0, 3.0, 30.0, PULSE_LONG, 103);
    ppi_post_antenna(10.0, 0.1, 1);
    ppi_post_plot(&a);
    ppi_post_antenna(100.0, 2.6, 1);
    ppi_post_plot(&c);
    ppi_post_antenna(210.0, 5.6, 1);
    ppi_post_plot(&b);
    ppi_post_antenna(215.0, 5.8, 1);
    w.refresh();
    check(w.plotCount() == 3, "three plots on the display");

    /* data tips */
    const QPointF ga = ground(a), gc = ground(c);
    QPointF sa = w.screenOf(ga.x(), ga.y());
    QString tip = w.tipText(sa + QPointF(3, 2));
    check(tip.contains("short pulse") && tip.contains("10.000 km") && tip.contains("12.50 m/s") &&
          tip.contains("no DBF") && tip.contains("Pulse number 101"), "hover data tip shows the plot values");
    check(w.tipText(sa + QPointF(30, 30)).isEmpty(), "no data tip away from the plots");
    send_mouse(w, QEvent::MouseMove, sa, Qt::NoButton);
    check(w.cursorText().contains("km") && w.cursorText().contains(QChar(0x00B0)), "cursor readout (range, azimuth)");

    /* wheel zoom around the cursor */
    const QPointF sc = w.screenOf(gc.x(), gc.y());
    send_wheel(w, sc, 120);
    const QPointF sc2 = w.screenOf(gc.x(), gc.y());
    check(std::fabs(w.zoom() - 1.25) < 1e-9 && std::hypot(sc2.x() - sc.x(), sc2.y() - sc.y()) < 0.5,
          "wheel zooms x1.25, the point under the cursor stays");
    for (int i = 0; i < 8; ++i) send_wheel(w, sc, 120);
    check(w.tipText(w.screenOf(gc.x(), gc.y())).contains("30.000 km"), "data tip works in the zoomed view");

    /* drag = pan */
    const QPointF v0 = w.viewCentre_m();
    const double k = w.screenOf(1000.0, 0.0).x() - w.screenOf(0.0, 0.0).x();   /* pixels per km */
    send_mouse(w, QEvent::MouseButtonPress, QPointF(400, 400), Qt::LeftButton);
    send_mouse(w, QEvent::MouseMove, QPointF(460, 400), Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, QPointF(460, 400), Qt::LeftButton);
    const QPointF v1 = w.viewCentre_m();
    check(std::fabs((v0.x() - v1.x()) - 60.0 / k * 1000.0) < 1.0 && std::fabs(v0.y() - v1.y()) < 1.0 &&
          w.pinnedCount() == 0, "drag pans the view (no pin)");

    /* click on a plot pins its data tip, a second click removes it */
    const QPointF sc3 = w.screenOf(gc.x(), gc.y());
    send_mouse(w, QEvent::MouseButtonPress, sc3, Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, sc3, Qt::LeftButton);
    check(w.pinnedCount() == 1, "click pins the data tip");
    if (!folder.empty()) w.saveSnapshot(QString::fromStdString(folder + "/ppi_zoom_pinned.png"));
    send_mouse(w, QEvent::MouseButtonPress, sc3, Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, sc3, Qt::LeftButton);
    check(w.pinnedCount() == 0, "second click removes it");

    /* double click and keys */
    send_mouse(w, QEvent::MouseButtonDblClick, QPointF(400, 400), Qt::LeftButton);
    check(w.zoom() == 1.0 && w.viewCentre_m() == QPointF(0, 0), "double click resets the view");
    send_key(w, Qt::Key_Plus, "+");
    check(std::fabs(w.zoom() - 1.5) < 1e-9, "key + zooms in");
    send_key(w, Qt::Key_Left, "");
    check(w.viewCentre_m().x() < 0.0, "arrow key pans");
    send_key(w, Qt::Key_R, "r");
    check(w.zoom() == 1.0, "key R resets the view");
    send_wheel(w, QPointF(400, 400), -120);
    check(w.zoom() == 1.0, "no zoom out beyond the full view");

    /* filters */
    w.setShowShort(false);
    sa = w.screenOf(ga.x(), ga.y());
    check(w.tipText(sa).isEmpty() && w.visibleCount() == 2, "short-pulse plots hidden (no data tip)");
    w.setShowShort(true);
    w.setMinSnr(20.0);
    check(w.visibleCount() == 2, "minimum SNR 20 dB hides the 14 dB plot");
    w.setMinSnr(-1e9);

    /* persistence: the antenna passes plot a (30 deg) again in scan 2 */
    w.setHistory(History::All);
    ppi_post_antenna(359.0, 9.9, 1);
    ppi_post_antenna(40.0, 11.2, 2);
    w.refresh();
    check(w.plotCount() == 3, "'all' keeps the plots");
    w.setHistory(History::UntilSwept);
    check(w.plotCount() == 2, "default: a plot disappears when the antenna passes it again");

    /* pause */
    w.setPaused(true);
    const Detection d = make_plot(20e3, 45.0, 2.0, 1.0, 20.0, PULSE_SHORT, 104);
    ppi_post_antenna(50.0, 11.3, 2);
    ppi_post_plot(&d);
    w.refresh();
    check(w.plotCount() == 2, "paused display does not take new plots");
    w.setPaused(false);
    w.refresh();
    check(w.plotCount() == 3, "resumed display shows them");
    w.clearPlots();
    check(w.plotCount() == 0, "clear plots");

    if (!folder.empty()) w.saveSnapshot(QString::fromStdString(folder + "/ppi_cleared.png"));
    std::printf(failures ? "%d FAILURE(S)\n" : "PASS\n", failures);
    return failures ? 1 : 0;
}
