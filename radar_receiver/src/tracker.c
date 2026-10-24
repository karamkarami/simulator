/*
 * tracker.c - extended Kalman filter tracking of the plots (as track_update.m).
 *
 * State relative to the antenna (x East, y North, z up): 2-D [x y vx vy],
 * 3-D [x y z vx vy vz], nearly constant velocity with white-acceleration
 * process noise q*[dt^3/3 dt^2/2; dt^2/2 dt] per axis. Measurement of a plot:
 * 2-D [ground range, azimuth] (slant range when the elevation is unknown),
 * 3-D [range, azimuth, elevation]; h(x) is linearised at the predicted state.
 * Association: smallest d^2 + ln|S| among the tracks whose chi-square gate
 * holds the plot; otherwise a new tentative track (zero velocity, velocity std
 * sigmaV0). Confirmation after confirmHits updates; deletion after more than
 * maxMiss (tentative: maxMissTentative) missed scans, a scan being missed
 * when the antenna passed the track half a scan ago without an update.
 */
#include "tracker.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

void tracker_init(Tracker *t, const RadarConfig *cfg)
{
    memset(t, 0, sizeof *t);
    t->dim = cfg->trackDim == 3 ? 3 : 2;
    t->n = 2 * t->dim;
    for (int i = 0; i < 3; i++) {
        t->R[i] = cfg->trackSigma[i] * cfg->trackSigma[i];
        t->q[i] = cfg->trackQ[i];
    }
    t->sigmaV0 = cfg->trackSigmaV0;
    t->gate = cfg->trackGate;
    t->scanTime_s = cfg->scanTime_s;
    t->height_m = cfg->height_m;
    t->confirmHits = cfg->trackConfirmHits;
    t->maxMissTentative = cfg->trackMaxMissTentative;
    t->maxMiss = cfg->trackMaxMiss;
    t->nextId = 1;
}

void tracker_free(Tracker *t)
{
    free(t->tracks);
    memset(t, 0, sizeof *t);
}

/* ---- small dense matrices (row-major) ----------------------------------- */

/* C = A (r x k) * B (k x c) */
static void mul(const double *A, const double *B, double *C, int r, int k, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            double s = 0.0;
            for (int m = 0; m < k; m++) s += A[i * k + m] * B[m * c + j];
            C[i * c + j] = s;
        }
}

/* C = A (r x k) * B' (B: c x k) */
static void mul_bt(const double *A, const double *B, double *C, int r, int k, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            double s = 0.0;
            for (int m = 0; m < k; m++) s += A[i * k + m] * B[j * k + m];
            C[i * c + j] = s;
        }
}

/* Inverse and determinant of a symmetric 2x2 or 3x3 matrix. Returns 0 if singular. */
static int inv_small(const double *S, int m, double *Si, double *det)
{
    if (m == 2) {
        *det = S[0] * S[3] - S[1] * S[2];
        if (!(*det > 0.0)) return 0;
        Si[0] = S[3] / *det;  Si[1] = -S[1] / *det;
        Si[2] = -S[2] / *det; Si[3] = S[0] / *det;
        return 1;
    }
    const double a = S[0], b = S[1], c = S[2], d = S[3], e = S[4], f = S[5], g = S[6], h = S[7], k = S[8];
    const double A = e * k - f * h, B = -(d * k - f * g), C = d * h - e * g;
    *det = a * A + b * B + c * C;
    if (!(*det > 0.0)) return 0;
    Si[0] = A / *det;  Si[1] = -(b * k - c * h) / *det; Si[2] = (b * f - c * e) / *det;
    Si[3] = B / *det;  Si[4] = (a * k - c * g) / *det;  Si[5] = -(a * f - c * d) / *det;
    Si[6] = C / *det;  Si[7] = -(a * h - b * g) / *det; Si[8] = (a * e - b * d) / *det;
    return 1;
}

/* ---- EKF ------------------------------------------------------------------ */

static void predict(const Tracker *t, const TrackState *s, double dt, double *x, double *P)
{
    const int d = t->dim, n = t->n;
    double F[36], FP[36];
    memset(F, 0, sizeof F);
    for (int i = 0; i < n; i++) F[i * n + i] = 1.0;
    for (int i = 0; i < d; i++) F[i * n + d + i] = dt;
    mul(F, s->x, x, n, n, 1);
    mul(F, s->P, FP, n, n, n);
    mul_bt(FP, F, P, n, n, n);
    const double a = fabs(dt);
    for (int i = 0; i < d; i++) {
        P[i * n + i]             += t->q[i] * a * a * a / 3.0;
        P[i * n + d + i]         += t->q[i] * a * a / 2.0;
        P[(d + i) * n + i]       += t->q[i] * a * a / 2.0;
        P[(d + i) * n + d + i]   += t->q[i] * a;
    }
}

/* Innovation z - h(x) (azimuth wrapped) and Jacobian H (dim x n). */
static void innovation(const Tracker *t, const double *z, const double *x, double *nu, double *H)
{
    const int n = t->n;
    const double px = x[0], py = x[1], rh2 = px * px + py * py, rh = sqrt(rh2);
    memset(H, 0, sizeof(double) * t->dim * n);
    if (t->dim == 3) {
        const double pz = x[2], R2 = rh2 + pz * pz, R = sqrt(R2);
        nu[0] = z[0] - R;
        nu[1] = z[1] - atan2(px, py);
        nu[2] = z[2] - atan2(pz, rh);
        H[0] = px / R;                    H[1] = py / R;                    H[2] = pz / R;
        H[n + 0] = py / rh2;              H[n + 1] = -px / rh2;
        H[2 * n + 0] = -px * pz / (R2 * rh); H[2 * n + 1] = -py * pz / (R2 * rh); H[2 * n + 2] = rh / R2;
    } else {
        nu[0] = z[0] - rh;
        nu[1] = z[1] - atan2(px, py);
        H[0] = px / rh;      H[1] = py / rh;
        H[n + 0] = py / rh2; H[n + 1] = -px / rh2;
    }
    nu[1] = fmod(nu[1] + PI, 2.0 * PI);                   /* MATLAB mod: result in [0, 2*pi) */
    if (nu[1] < 0.0) nu[1] += 2.0 * PI;
    nu[1] -= PI;
}

int tracker_plot(Tracker *t, const Detection *plot)
{
    const int d = t->dim, n = t->n;
    const double el = plot->el_deg * PI / 180.0, az = plot->az_deg * PI / 180.0;
    double z[3] = { 0.0, 0.0, 0.0 };
    if (d == 3) {
        if (!isfinite(el)) return -1;
        z[0] = plot->range_m;  z[1] = az;  z[2] = el;
    } else {
        z[0] = isfinite(el) ? plot->range_m * cos(el) : plot->range_m;
        z[1] = az;
    }
    const double tp = plot->time_s;

    int best = -1;
    double bestCost = INFINITY, xb[6], Pb[36], nub[3], Hb[18], Sib[9];
    for (int j = 0; j < t->nTracks; j++) {
        double x[6] = { 0.0 }, P[36], nu[3], H[18], HP[18], S[9], Si[9], det;
        predict(t, &t->tracks[j], tp - t->tracks[j].t, x, P);
        innovation(t, z, x, nu, H);
        mul(H, P, HP, d, n, n);
        mul_bt(HP, H, S, d, n, d);
        for (int i = 0; i < d; i++) S[i * d + i] += t->R[i];
        if (!inv_small(S, d, Si, &det)) continue;
        double Sinu[3], d2 = 0.0;
        mul(Si, nu, Sinu, d, d, 1);
        for (int i = 0; i < d; i++) d2 += nu[i] * Sinu[i];
        if (d2 > t->gate) continue;
        const double cost = d2 + log(det);
        if (cost < bestCost) {
            bestCost = cost;
            best = j;
            memcpy(xb, x, sizeof x);  memcpy(Pb, P, sizeof P);  memcpy(nub, nu, sizeof nu);
            memcpy(Hb, H, sizeof H);  memcpy(Sib, Si, sizeof Si);
        }
    }
    if (best >= 0) {                                       /* EKF update, Joseph form */
        double PHt[18], K[18], IKH[36], T1[36], P1[36], KR[18], P2[36];
        TrackState *s = &t->tracks[best];
        mul_bt(Pb, Hb, PHt, n, n, d);                     /* P H' (n x d) */
        mul(PHt, Sib, K, n, d, d);                         /* K = P H' S^-1 */
        for (int i = 0; i < n; i++) {
            double dx = 0.0;
            for (int m = 0; m < d; m++) dx += K[i * d + m] * nub[m];
            s->x[i] = xb[i] + dx;
        }
        mul(K, Hb, IKH, n, d, n);
        for (int i = 0; i < n * n; i++) IKH[i] = -IKH[i];
        for (int i = 0; i < n; i++) IKH[i * n + i] += 1.0;
        mul(IKH, Pb, T1, n, n, n);
        mul_bt(T1, IKH, P1, n, n, n);
        for (int i = 0; i < n; i++)
            for (int m = 0; m < d; m++) KR[i * d + m] = K[i * d + m] * t->R[m];
        mul_bt(KR, K, P2, n, d, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) s->P[i * n + j] = P1[i * n + j] + P2[i * n + j];
        for (int i = 0; i < n; i++)                        /* symmetrise */
            for (int j = i + 1; j < n; j++) {
                const double v = 0.5 * (s->P[i * n + j] + s->P[j * n + i]);
                s->P[i * n + j] = s->P[j * n + i] = v;
            }
        s->t = tp;
        s->hits++;
        s->misses = 0;
        s->confirmed = s->confirmed || s->hits >= t->confirmHits;
        return best;
    }

    /* new tentative track at the measured position */
    if (t->nTracks == t->capacity) {
        const int cap = t->capacity ? 2 * t->capacity : 32;
        TrackState *p = realloc(t->tracks, sizeof(TrackState) * cap);
        if (!p) return -1;
        t->tracks = p;
        t->capacity = cap;
    }
    TrackState *s = &t->tracks[t->nTracks];
    memset(s, 0, sizeof *s);
    double J[9], JR[9], Ppos[9];
    if (d == 3) {
        const double R = z[0], ce = cos(el), se = sin(el), sa = sin(az), ca = cos(az);
        s->x[0] = R * ce * sa;  s->x[1] = R * ce * ca;  s->x[2] = R * se;
        const double Jm[9] = { ce * sa, R * ce * ca, -R * se * sa,
                               ce * ca, -R * ce * sa, -R * se * ca,
                               se, 0.0, R * ce };
        memcpy(J, Jm, sizeof Jm);
    } else {
        const double rho = z[0], sa = sin(az), ca = cos(az);
        s->x[0] = rho * sa;  s->x[1] = rho * ca;
        const double Jm[4] = { sa, rho * ca, ca, -rho * sa };
        memcpy(J, Jm, sizeof Jm);
    }
    for (int i = 0; i < d; i++)
        for (int m = 0; m < d; m++) JR[i * d + m] = J[i * d + m] * t->R[m];
    mul_bt(JR, J, Ppos, d, d, d);                          /* J R J' */
    for (int i = 0; i < d; i++)
        for (int m = 0; m < d; m++) s->P[i * n + m] = Ppos[i * d + m];
    for (int i = d; i < n; i++) s->P[i * n + i] = t->sigmaV0 * t->sigmaV0;
    s->id = t->nextId++;
    s->t = tp;
    s->hits = 1;
    s->confirmed = t->confirmHits <= 1;
    return t->nTracks++;
}

void tracker_report(const Tracker *t, int index, Track *out)
{
    const TrackState *s = &t->tracks[index];
    const int d = t->dim;
    memset(out, 0, sizeof *out);
    out->id = s->id;
    out->confirmed = s->confirmed;
    out->dim = d;
    out->time_s = s->t;
    for (int i = 0; i < d; i++) {
        out->pos_m[i] = s->x[i];
        out->vel_mps[i] = s->x[d + i];
    }
    out->alt_m = d == 3 ? s->x[2] + t->height_m : NAN;
    out->hits = s->hits;
    out->misses = s->misses;
}

void tracker_maintain(Tracker *t, double t_s, void (*callback)(void *user, const Track *track), void *user)
{
    int keep = 0;
    for (int j = 0; j < t->nTracks; j++) {
        TrackState *s = &t->tracks[j];
        const double m = floor((t_s - s->t) / t->scanTime_s - 0.5);
        s->misses = m > 0.0 ? (int)m : 0;
        const int limit = s->confirmed ? t->maxMiss : t->maxMissTentative;
        if (s->misses > limit) {
            if (callback) {
                Track r;
                tracker_report(t, j, &r);
                r.deleted = 1;
                callback(user, &r);
            }
            continue;
        }
        if (keep != j) t->tracks[keep] = *s;
        keep++;
    }
    t->nTracks = keep;
}
