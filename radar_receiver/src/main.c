/*
 * main.c - radar receiver: DDC output over UDP -> processing -> plots.
 *
 *   radar_receiver [options]
 *     --config FILE          receiver configuration from export_receiver_config.m
 *                            (default config/radar_config.bin)
 *     --port N               UDP port of the DDC stream (default 50000)
 *     --record FILE          also store the received datagrams in FILE
 *     --replay FILE          read the datagrams from FILE instead of the network
 *     --ppi                  show the PPI (Qt build only)
 *     --ppi-snapshot FILE    save the PPI as PNG at the end and close it
 *     --serial DEVICE        send the CFAR detections to a serial port (e.g. /dev/ttyUSB0, COM3)
 *     --baud N               serial baud rate (default 921600)
 *     --serial-content X     'detections' (default) or 'plots'
 *     --threads N            processing threads (default: all cores but one)
 *     --update-every N       run Doppler/CFAR every N pulses (default from the configuration)
 *     --plots-csv FILE       plot list (default plots.csv, '' = none)
 *     --detections-csv FILE  list of all pruned CFAR detections (default: none)
 *     --idle-timeout S       stop after S seconds without data (network mode, default 0 = never)
 *     --quiet                do not print every plot
 */
#include "radar_rx.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *configPath, *recordPath, *replayPath, *snapshotPath;
    const char *serialDevice, *serialContent, *plotsCsv, *detectionsCsv;
    int port, baud, threads, updateEvery, ppi, quiet;
    double idleTimeout_s;
} Options;

typedef struct {
    Options opt;
    RadarConfig cfg;
    SerialOut *serial;
    FILE *plotsFile, *detFile;
    uint64_t nPlots, nDetections;
    int exitCode;
} App;

static volatile int g_stop = 0;

static void on_signal(int sig)
{
    (void)sig;
    g_stop = 1;
}

static void usage(void)
{
    puts("usage: radar_receiver [--config FILE] [--port N] [--record FILE] [--replay FILE]\n"
         "                      [--ppi] [--ppi-snapshot FILE] [--serial DEVICE] [--baud N]\n"
         "                      [--serial-content detections|plots] [--threads N] [--update-every N]\n"
         "                      [--plots-csv FILE] [--detections-csv FILE] [--idle-timeout S] [--quiet]");
}

static int parse_options(int argc, char **argv, Options *o)
{
    memset(o, 0, sizeof(*o));
    o->configPath = "config/radar_config.bin";
    o->port = 50000;
    o->baud = 921600;
    o->serialContent = "detections";
    o->plotsCsv = "plots.csv";
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
#define TAKE(field) do { if (!v) { fprintf(stderr, "%s needs a value\n", a); return -1; } field; i++; } while (0)
        if      (!strcmp(a, "--config"))         TAKE(o->configPath = v);
        else if (!strcmp(a, "--port"))           TAKE(o->port = atoi(v));
        else if (!strcmp(a, "--record"))         TAKE(o->recordPath = v);
        else if (!strcmp(a, "--replay"))         TAKE(o->replayPath = v);
        else if (!strcmp(a, "--ppi"))            o->ppi = 1;
        else if (!strcmp(a, "--ppi-snapshot"))   TAKE((o->snapshotPath = v, o->ppi = 1));
        else if (!strcmp(a, "--serial"))         TAKE(o->serialDevice = v);
        else if (!strcmp(a, "--baud"))           TAKE(o->baud = atoi(v));
        else if (!strcmp(a, "--serial-content")) TAKE(o->serialContent = v);
        else if (!strcmp(a, "--threads"))        TAKE(o->threads = atoi(v));
        else if (!strcmp(a, "--update-every"))   TAKE(o->updateEvery = atoi(v));
        else if (!strcmp(a, "--plots-csv"))      TAKE(o->plotsCsv = v);
        else if (!strcmp(a, "--detections-csv")) TAKE(o->detectionsCsv = v);
        else if (!strcmp(a, "--idle-timeout"))   TAKE(o->idleTimeout_s = atof(v));
        else if (!strcmp(a, "--quiet"))          o->quiet = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); return 1; }
        else { fprintf(stderr, "unknown option %s\n", a); usage(); return -1; }
#undef TAKE
    }
    if (strcmp(o->serialContent, "detections") && strcmp(o->serialContent, "plots")) {
        fprintf(stderr, "--serial-content must be 'detections' or 'plots'\n");
        return -1;
    }
    return 0;
}

/* ------------------------------------------------------------------------ */
/* Callbacks of the processor                                               */
/* ------------------------------------------------------------------------ */

static const char *pulse_name(int pulse) { return pulse == PULSE_SHORT ? "short" : "long"; }

static void on_detections(void *user, int pulse, const Detection *dets, int n)
{
    App *app = user;
    app->nDetections += (uint64_t)n;
    /* called on every CPI update, also without detections: no frame is sent
     * then, but the bytes still queued are written to the port */
    if (app->serial) {
        if (n > 0 && !strcmp(app->opt.serialContent, "detections"))
            serial_send_detections(app->serial, pulse, dets[0].pulseNumber, dets[0].az_deg, dets, n);
        else
            serial_send_detections(app->serial, pulse, 0, 0.0, NULL, 0);
    }
    if (app->detFile) {
        for (int i = 0; i < n; i++) {
            const Detection *d = &dets[i];
            fprintf(app->detFile, "%u,%s,%.6f,%.3f,%.1f,%.3f,%.3f,%.2f,%d,%d,%d\n", d->pulseNumber, pulse_name(pulse),
                    d->time_s, d->az_deg, d->range_m, d->el_deg, d->vr_mps, d->snr_dB, d->beam, d->rIdx, d->dIdx);
        }
    }
}

static void on_plot(void *user, const Detection *p)
{
    App *app = user;
    app->nPlots++;
    if (!app->opt.quiet)
        printf("plot  t %8.3f s  az %7.2f deg  R %8.3f km  el %6.2f deg  vr %7.2f m/s  SNR %5.1f dB  %-5s  hits %d\n",
               p->time_s, p->az_deg, p->range_m / 1e3, p->el_deg, p->vr_mps, p->snr_dB, pulse_name(p->pulse), p->hits);
    if (app->plotsFile) {
        fprintf(app->plotsFile, "%.6f,%.4f,%.2f,%.4f,%.4f,%.2f,%s,%d,%d,%d,%d,%u\n", p->time_s, p->az_deg, p->range_m,
                p->el_deg, p->vr_mps, p->snr_dB, pulse_name(p->pulse), p->hits, p->beam, p->rIdx, p->dIdx, p->pulseNumber);
        fflush(app->plotsFile);
    }
    if (app->serial && !strcmp(app->opt.serialContent, "plots")) serial_send_plot(app->serial, p);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) ppi_post_plot(p);
#endif
}

static void on_antenna(void *user, double az_deg, double time_s, uint32_t scan)
{
    App *app = user;
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) ppi_post_antenna(az_deg, time_s, scan);
#else
    (void)app; (void)az_deg; (void)time_s; (void)scan;
#endif
}

/* ------------------------------------------------------------------------ */
/* Worker: input -> processing                                              */
/* ------------------------------------------------------------------------ */

static void print_timing(const Processor *proc, const UdpInput *in, double wall_s)
{
    ProcessorTiming t;
    UdpStats u;
    processor_timing(proc, &t);
    udp_input_stats(in, &u);
    const double n = t.pris ? (double)t.pris : 1.0;
    fprintf(stderr, "\nPRIs processed: %llu (%llu CPI updates) in %.1f s; datagrams %llu, foreign %llu, "
            "incomplete PRIs %llu\n", (unsigned long long)t.pris, (unsigned long long)t.updates, wall_s,
            (unsigned long long)u.datagrams, (unsigned long long)u.foreign, (unsigned long long)u.prisIncomplete);
    fprintf(stderr, "Processing time per PRI: %.3f ms (DBF %.3f, MF %.3f, MTI %.3f, Doppler %.3f, "
            "CFAR+pruning %.3f, extraction %.3f)\n", t.total_ms / n, t.dbf_ms / n, t.mf_ms / n,
            t.mti_ms / n, t.doppler_ms / n, t.cfar_ms / n, t.extract_ms / n);
}

static int worker(void *arg)
{
    App *app = arg;
    char err[256] = "";
    UdpInput *in = udp_input_open(app->opt.port, app->opt.recordPath, app->opt.replayPath, err, sizeof err);
    if (!in) {
        fprintf(stderr, "input: %s\n", err);
        app->exitCode = 1;
        return 1;
    }
    Processor *proc = processor_create(&app->cfg, app->opt.threads);
    if (!proc) {
        fprintf(stderr, "cannot create the processor\n");
        udp_input_close(in);
        app->exitCode = 1;
        return 1;
    }
    ProcessorCallbacks cb = { on_detections, on_plot, on_antenna, app };
    PriData pri;
    memset(&pri, 0, sizeof pri);
    if (!app->opt.replayPath)
        fprintf(stderr, "Waiting for the DDC stream on UDP port %d ...\n", app->opt.port);
    const double t0 = now_seconds();
    double lastData = t0;
    while (!g_stop) {
        int r = udp_input_next(in, &pri, 200);
        if (r < 0) break;                                      /* end of the replay file */
        if (r == 0) {
            if (app->opt.idleTimeout_s > 0 && now_seconds() - lastData > app->opt.idleTimeout_s) break;
            continue;
        }
        lastData = now_seconds();
        processor_process(proc, &pri, &cb);
    }
    processor_flush(proc, &cb);
    print_timing(proc, in, now_seconds() - t0);
    fprintf(stderr, "Plots: %llu, CFAR detections after pruning: %llu\n",
            (unsigned long long)app->nPlots, (unsigned long long)app->nDetections);
    processor_destroy(proc);
    pri_data_free(&pri);
    udp_input_close(in);
    return 0;
}

int main(int argc, char **argv)
{
    App app;
    memset(&app, 0, sizeof app);
    int r = parse_options(argc, argv, &app.opt);
    if (r != 0) return r < 0 ? 2 : 0;

    char err[256] = "";
    if (radar_config_load(app.opt.configPath, &app.cfg, err, sizeof err) != 0) {
        fprintf(stderr, "configuration: %s\n", err);
        return 1;
    }
    if (app.opt.updateEvery > 0 && app.opt.updateEvery != app.cfg.updateEvery) {
        /* the plot extractor needs the same azimuth extent with the new update interval */
        app.cfg.updateEvery = app.opt.updateEvery;
        const double step_deg = app.cfg.rate_dps * app.cfg.meanPri_s * app.cfg.updateEvery;
        const int minHits = (int)ceil(app.cfg.minExtent_deg / step_deg);     /* as run_radar_sim.m */
        app.cfg.minHits = minHits > 1 ? minHits : 1;
    }
    fprintf(stderr, "Configuration: %d element(s), %d beam(s), DDC %d samples, short %d / long %d range cells, "
            "CPI %d pulses, FFT %d x %d row(s), update every %d pulse(s), plots need %d update(s)\n", app.cfg.nCh,
            app.cfg.nBeams, app.cfg.nDdc, app.cfg.pulse[0].nRange, app.cfg.pulse[1].nRange, app.cfg.cpiLength,
            app.cfg.nfft, app.cfg.nRows, app.cfg.updateEvery, app.cfg.minHits);

    if (app.opt.plotsCsv && app.opt.plotsCsv[0]) {
        app.plotsFile = fopen(app.opt.plotsCsv, "w");
        if (app.plotsFile)
            fprintf(app.plotsFile, "time_s,az_deg,range_m,el_deg,vr_mps,snr_dB,pulse,hits,beam,rangeCell,dopplerBin,pulseNumber\n");
    }
    if (app.opt.detectionsCsv && app.opt.detectionsCsv[0]) {
        app.detFile = fopen(app.opt.detectionsCsv, "w");
        if (app.detFile)
            fprintf(app.detFile, "pulseNumber,pulse,time_s,az_deg,range_m,el_deg,vr_mps,snr_dB,beam,rangeCell,dopplerBin\n");
    }
    if (app.opt.serialDevice) {
        app.serial = serial_open(app.opt.serialDevice, app.opt.baud, err, sizeof err);
        if (!app.serial) {
            fprintf(stderr, "serial: %s\n", err);
            return 1;
        }
    }
    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    if (app.opt.ppi) {
#ifdef RADAR_WITH_QT
        app.exitCode = ppi_run(argc, argv, app.cfg.maxRange_m, worker, &app, &g_stop, app.opt.snapshotPath);
#else
        fprintf(stderr, "this build has no Qt: --ppi is not available\n");
        app.exitCode = 1;
#endif
    } else {
        worker(&app);
    }

    if (app.serial) {
        uint64_t sent, dropped, bytes;
        serial_stats(app.serial, &sent, &dropped, &bytes);
        fprintf(stderr, "Serial: %llu frames sent, %llu dropped, %llu bytes\n", (unsigned long long)sent,
                (unsigned long long)dropped, (unsigned long long)bytes);
        serial_close(app.serial);
    }
    if (app.plotsFile) fclose(app.plotsFile);
    if (app.detFile) fclose(app.detFile);
    radar_config_free(&app.cfg);
    return app.exitCode;
}
