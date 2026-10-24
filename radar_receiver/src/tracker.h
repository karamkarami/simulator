/*
 * tracker.h - EKF tracker of the plots (internal).
 * Same algorithm as src/track_update.m of the MATLAB simulator.
 */
#ifndef TRACKER_H
#define TRACKER_H

#include "radar_rx.h"

typedef struct {
    uint32_t id;
    int      confirmed;
    double   x[6];          /* [x y vx vy] (2-D) or [x y z vx vy vz] (3-D)      */
    double   P[36];         /* covariance, row-major n x n                      */
    double   t;             /* time of the last update                          */
    int      hits, misses;
} TrackState;

typedef struct {
    int      dim, n;                     /* measurement dimension, state size 2*dim */
    double   R[3];                       /* measurement variances                   */
    double   q[3];
    double   sigmaV0, gate, scanTime_s, height_m;
    int      confirmHits, maxMissTentative, maxMiss;
    TrackState *tracks;
    int      nTracks, capacity;
    uint32_t nextId;
} Tracker;

void tracker_init(Tracker *t, const RadarConfig *cfg);
void tracker_free(Tracker *t);
/* Associates the plot with a track (EKF update) or starts a new tentative
 * track. Returns the index of that track, or -1 when the plot is not usable
 * (3-D without elevation). */
int  tracker_plot(Tracker *t, const Detection *plot);
/* Counts the missed scans at time t_s and deletes the tracks with too many;
 * every deleted track is reported through callback (deleted = 1). */
void tracker_maintain(Tracker *t, double t_s, void (*callback)(void *user, const Track *track), void *user);
void tracker_report(const Tracker *t, int index, Track *out);

#endif
