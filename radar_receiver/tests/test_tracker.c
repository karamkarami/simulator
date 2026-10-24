/*
 * test_tracker.c - the C EKF tracker (tracker.c) against the MATLAB tracker.
 *
 *   test_tracker CONFIG                      self test: two crossing targets, 3-D and 2-D
 *   test_tracker CONFIG PLOTS.csv LOG.csv DIM GATE
 *       PLOTS.csv: time_s,az_deg,range_m,el_deg of every plot in processing order
 *       LOG.csv:   the MATLAB track log (track_update.m): time,id,confirmed,x,y,z,vx,vy,vz,hits,plotIndex
 *       DIM, GATE: tracker dimension and chi-square gate (they replace those of CONFIG)
 *       Every plot is given to the tracker (with maintenance at the plot time,
 *       as tests/test_tracker.m does); every track update must match the log:
 *       same track number, state and number of updates, positions and
 *       velocities within 1e-6 relative.
 *
 * Prints PASS or FAIL and returns 0 on success.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "radar_rx.h"
#include "tracker.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
static int failures = 0;

#define CHECK(cond, ...)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL: " __VA_ARGS__);                                    \
            printf("\n");                                                    \
            failures++;                                                      \
        }                                                                    \
    } while (0)

static double gauss(void)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0), u2 = rand() / (RAND_MAX + 1.0);
    return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

/* Two targets crossing each other, 6 scans, plot errors as configured. */
static void self_test(RadarConfig *cfg, int dim)
{
    cfg->trackDim = dim;
    Tracker t;
    tracker_init(&t, cfg);
    const double p0[2][3] = { { 20e3, 5e3, 2000 }, { 20e3, -5e3, 3000 } };
    const double v[2][3] = { { -50, -150, 0 }, { -50, 150, 0 } };
    const double T = cfg->scanTime_s, h = cfg->height_m;
    srand(3);
    for (int scan = 0; scan < 6; scan++) {
        for (int i = 0; i < 2; i++) {
            const double tt = scan * T + i * 0.1;
            const double x = p0[i][0] + v[i][0] * tt, y = p0[i][1] + v[i][1] * tt, z = p0[i][2] - h;
            const double R = sqrt(x * x + y * y + z * z);
            Detection d;
            memset(&d, 0, sizeof d);
            d.time_s = tt;
            d.range_m = R + cfg->trackSigma[0] * gauss();
            d.az_deg = (atan2(x, y) + cfg->trackSigma[1] * gauss()) * 180.0 / PI;
            d.el_deg = dim == 3 ? (asin(z / R) + cfg->trackSigma[2] * gauss()) * 180.0 / PI : NAN;
            tracker_plot(&t, &d);
            tracker_maintain(&t, tt, NULL, NULL);
        }
    }
    int confirmed = 0;
    for (int j = 0; j < t.nTracks; j++) {
        Track r;
        tracker_report(&t, j, &r);
        if (!r.confirmed) continue;
        confirmed++;
        const double tt = r.time_s;
        double dist[2];
        for (int k = 0; k < 2; k++)
            dist[k] = hypot(r.pos_m[0] - (p0[k][0] + v[k][0] * tt), r.pos_m[1] - (p0[k][1] + v[k][1] * tt));
        const int i = dist[0] <= dist[1] ? 0 : 1;          /* nearest target (they cross at t = 33 s) */
        const double ex = r.pos_m[0] - (p0[i][0] + v[i][0] * tt), ey = r.pos_m[1] - (p0[i][1] + v[i][1] * tt);
        const double ev = hypot(r.vel_mps[0] - v[i][0], r.vel_mps[1] - v[i][1]);
        printf("  %d-D: track %u, %d updates, position error %.0f m, velocity error %.1f m/s", dim, r.id, r.hits,
               hypot(ex, ey), ev);
        if (dim == 3) printf(", height %.0f m (true %.0f)", r.alt_m, p0[i][2]);
        printf("\n");
        CHECK(r.hits == 6, "%d-D: track %u has %d updates (expected 6)", dim, r.id, r.hits);
        if (dim == 3) CHECK(hypot(ex, ey) < 300 && ev < 30, "%d-D: track %u too far from its target", dim, r.id);
        if (dim == 3) CHECK(fabs(r.alt_m - p0[i][2]) < 400, "%d-D: height of track %u", dim, r.id);
    }
    CHECK(confirmed == 2, "%d-D: %d confirmed tracks (expected 2)", dim, confirmed);
    tracker_free(&t);
}

static int compare_with_matlab(RadarConfig *cfg, const char *plotsPath, const char *logPath, int dim, double gate)
{
    FILE *fp = fopen(plotsPath, "r"), *fl = fopen(logPath, "r");
    if (!fp || !fl) {
        printf("FAIL: cannot open %s or %s\n", plotsPath, logPath);
        if (fp) fclose(fp);
        if (fl) fclose(fl);
        return 0;
    }
    cfg->trackDim = dim;
    cfg->trackGate = gate;
    Tracker t;
    tracker_init(&t, cfg);
    double ref[11];
    int haveRef = fscanf(fl, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &ref[0], &ref[1], &ref[2], &ref[3], &ref[4],
                         &ref[5], &ref[6], &ref[7], &ref[8], &ref[9], &ref[10]) == 11;
    int plotIndex = 0, updates = 0, bad = 0;
    double maxRel = 0.0;
    double tp, az, R, el;
    while (fscanf(fp, "%lf,%lf,%lf,%lf", &tp, &az, &R, &el) == 4) {
        plotIndex++;
        Detection d;
        memset(&d, 0, sizeof d);
        d.time_s = tp;  d.az_deg = az;  d.range_m = R;  d.el_deg = el;
        const int i = tracker_plot(&t, &d);
        if (i >= 0) {
            Track r;
            tracker_report(&t, i, &r);
            updates++;
            if (!haveRef || (int)ref[10] != plotIndex) {
                bad++;
                printf("FAIL: plot %d updated track %u, MATLAB did not\n", plotIndex, r.id);
            } else {
                const double c[9] = { r.time_s, (double)r.id, (double)r.confirmed, r.pos_m[0], r.pos_m[1], r.pos_m[2],
                                      r.vel_mps[0], r.vel_mps[1], r.vel_mps[2] };
                int same = c[1] == ref[1] && c[2] == ref[2] && r.hits == (int)ref[9];
                for (int k = 3; k < 9; k++) {
                    const double rel = fabs(c[k] - ref[k]) / (fabs(ref[k]) + 1.0);
                    if (rel > maxRel) maxRel = rel;
                    same &= rel < 1e-6;
                }
                if (!same) {
                    bad++;
                    if (bad < 5) printf("FAIL: plot %d: C track %u (conf %d, hits %d), MATLAB track %.0f (conf %.0f, hits %.0f)\n",
                                        plotIndex, r.id, r.confirmed, r.hits, ref[1], ref[2], ref[9]);
                }
                haveRef = fscanf(fl, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &ref[0], &ref[1], &ref[2], &ref[3],
                                 &ref[4], &ref[5], &ref[6], &ref[7], &ref[8], &ref[9], &ref[10]) == 11;
            }
        }
        tracker_maintain(&t, tp, NULL, NULL);
    }
    if (haveRef) { bad++; printf("FAIL: MATLAB has more track updates than C\n"); }
    printf("  %d-D: %d plots, %d track updates, %d different from MATLAB, max relative state difference %.1e\n",
           dim, plotIndex, updates, bad, maxRel);
    fclose(fp);
    fclose(fl);
    tracker_free(&t);
    failures += bad > 0;
    return bad == 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: test_tracker CONFIG [PLOTS.csv LOG.csv DIM GATE]\n");
        return 2;
    }
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(argv[1], &cfg, err, sizeof err) != 0) {
        printf("FAIL: configuration %s: %s\n", argv[1], err);
        return 1;
    }
    if (!cfg.trackEnable) printf("(tracking is off in the configuration; its parameters are used anyway)\n");
    if (argc >= 6) {
        compare_with_matlab(&cfg, argv[2], argv[3], atoi(argv[4]), atof(argv[5]));
    } else {
        self_test(&cfg, 3);
        self_test(&cfg, 2);
    }
    radar_config_free(&cfg);
    printf(failures ? "FAIL\n" : "PASS\n");
    return failures ? 1 : 0;
}
