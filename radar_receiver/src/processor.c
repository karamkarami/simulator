/*
 * processor.c - receiver processing chain of one PRI, for the short and the
 * long sub-pulse separately (same algorithms as the MATLAB simulator):
 *
 *   DDC output (int16 I/Q) -> DBF -> matched filter -> CPI buffer -> MTI
 *   -> Doppler (window, FFT, non-coherent integration of rows, noise-gain
 *   equalisation) -> 2-D CA-CFAR -> pruning -> detection reports -> plot
 *   extraction.
 *
 * Speed: long matched filters run as fast convolution (FFT, identical result
 * to the valid convolution); the MTI output of every pulse is computed once
 * and kept in a ring buffer; the Doppler FFTs of all beam/range cells run
 * together (cells are the fast index, so the loops vectorise) in cache-sized
 * blocks; the CFAR works beam by beam with running Doppler sums so its data
 * stays in the cache; DBF, matched filter, Doppler and CFAR are split over
 * threads with OpenMP. Work buffers are allocated once.
 *
 * Array layouts (float, real and imaginary parts in separate arrays):
 *   cell             = beam * nRange + rangeCell
 *   CPI / MTI buffer = [slot][cell]
 *   power map P      = [dopplerBin][cell]   (Doppler after fftshift)
 */
#include "radar_rx.h"
#include "extractor.h"
#include "fft.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.14159265358979323846
#define CELL_BLOCK 64           /* cells per Doppler FFT block (fits the L1/L2 cache) */
#define MF_LANES   24           /* beams per fast-convolution chunk (3 AVX registers) */

typedef struct {
    int nR, nCells, nSlots, nMtiOut;
    float *segRe, *segIm;       /* DDC samples of this sub-pulse [ch][n]          */
    float *beamRe, *beamIm;     /* beamformer output [beam][n]                    */
    float *cpiRe, *cpiIm;       /* matched filter outputs [slot][cell]            */
    float *mtiRe, *mtiIm;       /* MTI outputs [slot][cell]                       */
    int    last, count, nPushed;            /* CPI buffer state                   */
    int    mtiLast, mtiCount;               /* MTI ring state                     */
    uint16_t *staggerIdx;       /* stagger position of every buffered pulse        */
    double   *az, *time;
    uint32_t newestPulse;
    float *P;                   /* power map [bin][cell]                          */
    unsigned char *det;         /* CFAR detections [bin][cell]                    */
    Detection *reports;         /* pruned detections of the latest CPI           */
    int nReports, capReports;
    Extractor ext;
    /* matched filter by fast convolution (used when cheaper than direct) */
    int    mfFft, mfN, mfChunks;
    FftPlan mfPlan;
    float *mfHre, *mfHim;       /* FFT of the zero-padded taps, divided by mfN    */
    float *mfWork;              /* per chunk: Xre, Xim, Yre, Yim [mfN][MF_LANES]  */
} PulseState;

typedef struct {                /* work buffers of one thread                      */
    float  *xr, *xi;            /* Doppler: one cell block [nfft][CELL_BLOCK]     */
    float  *sumO, *sumG;        /* CFAR: range box sums of one beam [nfft][nR]    */
    double *prefix;             /* CFAR: prefix sums of four range rows           */
    double *accO, *accG;        /* CFAR: Doppler running sums [nR]                */
} Scratch;

struct Processor {
    const RadarConfig *cfg;
    FftPlan fft;
    PulseState ps[N_PULSES];
    int nThreads;
    Scratch *scratch;           /* [nThreads] */
    uint32_t lastPulse;
    int havePulse;
    ProcessorTiming timing;
};

/* ------------------------------------------------------------------------ */

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) abort();
    return p;
}

static void *xrealloc(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) abort();
    return p;
}

static void reset_buffers(PulseState *s)
{
    s->last = -1;
    s->count = 0;
    s->nPushed = 0;
    s->mtiLast = -1;
    s->mtiCount = 0;
}

static int thread_num(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* Fast convolution of the matched filter: FFT length and the FFT of the taps.
 * It is used when its estimated time is shorter than that of the direct
 * convolution (which is split over the threads beam by beam; the FFTs of
 * MF_LANES beams run in one thread). Costs per operation measured on an
 * AVX2 Xeon: 0.3 ns per complex MAC (direct), 1.5 ns per FFT unit. */
static int mf_fft_init(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int forceDirect, int nThreads)
{
    int N = 1;
    while (N < pc->ddcCount) N *= 2;
    int log2N = 0;
    while ((1 << log2N) < N) log2N++;
    const int nChunks = (cfg->nBeams + MF_LANES - 1) / MF_LANES;
    const double direct = 0.3 * pc->nTaps * pc->nRange * cfg->nBeams / nThreads;
    const double fast = 1.5 * ((double)N * log2N + N) * nChunks * MF_LANES / (nThreads < nChunks ? nThreads : nChunks);
    s->mfFft = !forceDirect && fast < direct;
    if (!s->mfFft) return 0;
    s->mfN = N;
    s->mfChunks = nChunks;
    if (fft_plan_init(&s->mfPlan, N) != 0) return -1;
    s->mfHre = xcalloc(N, sizeof(float));
    s->mfHim = xcalloc(N, sizeof(float));
    for (int j = 0; j < pc->nTaps; j++) {                           /* bit-reversed input      */
        s->mfHre[s->mfPlan.bitrev[j]] = pc->tapsRe[j] / N;
        s->mfHim[s->mfPlan.bitrev[j]] = pc->tapsIm[j] / N;
    }
    fft_rows(&s->mfPlan, s->mfHre, s->mfHim, 1, 1);
    s->mfWork = xcalloc((size_t)s->mfChunks * 4 * N * MF_LANES, sizeof(float));
    return 0;
}

Processor *processor_create(const RadarConfig *cfg, int nThreads)
{
#ifdef _OPENMP
    /* default: all cores but one (left for the UDP input and the display) */
    if (nThreads <= 0) nThreads = omp_get_num_procs() > 1 ? omp_get_num_procs() - 1 : 1;
    omp_set_num_threads(nThreads);
    const int nThr = omp_get_max_threads();
#else
    (void)nThreads;
    const int nThr = 1;
#endif
    Processor *p = xcalloc(1, sizeof(Processor));
    p->cfg = cfg;
    if (fft_plan_init(&p->fft, cfg->nfft) != 0) { free(p); return NULL; }
    int maxNR = 1;
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        s->nR = pc->nRange;
        if (s->nR > maxNR) maxNR = s->nR;
        s->nCells = cfg->nBeams * s->nR;
        s->nSlots = cfg->cpiLength;
        s->nMtiOut = cfg->cpiLength - cfg->nMtiTaps + 1;
        size_t seg = (size_t)cfg->nCh * pc->ddcCount, beams = (size_t)cfg->nBeams * pc->ddcCount;
        size_t cpi = (size_t)s->nSlots * s->nCells, mti = (size_t)s->nMtiOut * s->nCells;
        size_t map = (size_t)cfg->nfft * s->nCells;
        s->segRe = xcalloc(seg, sizeof(float));   s->segIm = xcalloc(seg, sizeof(float));
        s->beamRe = xcalloc(beams, sizeof(float)); s->beamIm = xcalloc(beams, sizeof(float));
        s->cpiRe = xcalloc(cpi, sizeof(float));   s->cpiIm = xcalloc(cpi, sizeof(float));
        s->mtiRe = xcalloc(mti, sizeof(float));   s->mtiIm = xcalloc(mti, sizeof(float));
        s->staggerIdx = xcalloc(s->nSlots, sizeof(uint16_t));
        s->az = xcalloc(s->nSlots, sizeof(double));
        s->time = xcalloc(s->nSlots, sizeof(double));
        s->P = xcalloc(map, sizeof(float));
        s->det = xcalloc(map, 1);
        if (mf_fft_init(cfg, pc, s, cfg->mfDirect, nThr) != 0) { processor_destroy(p); return NULL; }
        extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
        reset_buffers(s);
    }
    p->nThreads = nThr;
    p->scratch = xcalloc(nThr, sizeof(Scratch));
    for (int t = 0; t < nThr; t++) {
        Scratch *w = &p->scratch[t];
        w->xr = xcalloc((size_t)cfg->nfft * CELL_BLOCK, sizeof(float));
        w->xi = xcalloc((size_t)cfg->nfft * CELL_BLOCK, sizeof(float));
        w->sumO = xcalloc((size_t)cfg->nfft * maxNR, sizeof(float));
        w->sumG = xcalloc((size_t)cfg->nfft * maxNR, sizeof(float));
        w->prefix = xcalloc(4 * ((size_t)maxNR + 1), sizeof(double));
        w->accO = xcalloc(maxNR, sizeof(double));
        w->accG = xcalloc(maxNR, sizeof(double));
    }
    return p;
}

void processor_destroy(Processor *p)
{
    if (!p) return;
    for (int i = 0; i < N_PULSES; i++) {
        PulseState *s = &p->ps[i];
        free(s->segRe); free(s->segIm); free(s->beamRe); free(s->beamIm);
        free(s->cpiRe); free(s->cpiIm); free(s->mtiRe); free(s->mtiIm);
        free(s->staggerIdx); free(s->az); free(s->time);
        free(s->P); free(s->det);
        free(s->reports);
        if (s->mfFft) fft_plan_free(&s->mfPlan);
        free(s->mfHre); free(s->mfHim); free(s->mfWork);
        extractor_free(&s->ext);
    }
    for (int t = 0; p->scratch && t < p->nThreads; t++) {
        Scratch *w = &p->scratch[t];
        free(w->xr); free(w->xi); free(w->sumO); free(w->sumG); free(w->prefix); free(w->accO); free(w->accG);
    }
    free(p->scratch);
    fft_plan_free(&p->fft);
    free(p);
}

void processor_timing(const Processor *p, ProcessorTiming *t)
{
    *t = p->timing;
}

/* ------------------------------------------------------------------------ */
/* Fast time: DBF and matched filter                                        */
/* ------------------------------------------------------------------------ */

/* DDC samples of the sub-pulse, as float [ch][n]. */
static void load_segment(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s)
{
    const int n = pc->ddcCount;
    for (int ch = 0; ch < cfg->nCh; ch++) {
        const int16_t *src = pri->iq + 2 * ((size_t)ch * pri->nSamples + pc->ddcFirst);
        float *re = s->segRe + (size_t)ch * n, *im = s->segIm + (size_t)ch * n;
        for (int k = 0; k < n; k++) {
            re[k] = src[2 * k];
            im[k] = src[2 * k + 1];
        }
    }
}

/* beam b = sum over channels of W[b][ch] * x[ch]  (or the single channel). */
static void beamform(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s)
{
    const int n = pc->ddcCount;
    if (!cfg->dbfEnabled) {
        memcpy(s->beamRe, s->segRe, sizeof(float) * n);
        memcpy(s->beamIm, s->segIm, sizeof(float) * n);
        return;
    }
    #pragma omp parallel for schedule(static)
    for (int b = 0; b < cfg->nBeams; b++) {
        float *restrict yr = s->beamRe + (size_t)b * n, *restrict yi = s->beamIm + (size_t)b * n;
        memset(yr, 0, sizeof(float) * n);
        memset(yi, 0, sizeof(float) * n);
        for (int ch = 0; ch < cfg->nCh; ch++) {
            const float wr = cfg->wRe[b * cfg->nCh + ch], wi = cfg->wIm[b * cfg->nCh + ch];
            const float *restrict xr = s->segRe + (size_t)ch * n, *restrict xi = s->segIm + (size_t)ch * n;
            for (int k = 0; k < n; k++) {
                yr[k] += wr * xr[k] - wi * xi[k];
                yi[k] += wr * xi[k] + wi * xr[k];
            }
        }
    }
}

/* The same valid convolution by FFT: the circular convolution of length
 * mfN >= ddcCount equals the linear one at the valid outputs r + nTaps - 1.
 * Inverse FFT = conj(FFT(conj(.))) (1/mfN is in the taps). The beams are
 * processed in chunks of MF_LANES that are the fast index of the FFT rows. */
static void matched_filter_fft(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot)
{
    const int N = s->mfN, n = pc->ddcCount, L = pc->nTaps, nR = s->nR, nB = cfg->nBeams;
    const int *rev = s->mfPlan.bitrev;
    #pragma omp parallel for schedule(static)
    for (int c = 0; c < s->mfChunks; c++) {
        const int b0 = c * MF_LANES, nb = (nB - b0 < MF_LANES) ? nB - b0 : MF_LANES;
        float *Xr = s->mfWork + (size_t)c * 4 * N * MF_LANES;
        float *Xi = Xr + (size_t)N * MF_LANES, *Yr = Xi + (size_t)N * MF_LANES, *Yi = Yr + (size_t)N * MF_LANES;
        for (int m = 0; m < N; m++) {                           /* beams -> bit-reversed rows */
            float *dr = Xr + (size_t)rev[m] * MF_LANES, *di = Xi + (size_t)rev[m] * MF_LANES;
            for (int j = 0; j < MF_LANES; j++) {
                const int ok = m < n && j < nb;
                dr[j] = ok ? s->beamRe[(size_t)(b0 + j) * n + m] : 0.0f;
                di[j] = ok ? s->beamIm[(size_t)(b0 + j) * n + m] : 0.0f;
            }
        }
        fft_rows(&s->mfPlan, Xr, Xi, MF_LANES, MF_LANES);
        for (int k = 0; k < N; k++) {                           /* conj(X * H), bit-reversed  */
            const float hr = s->mfHre[k], hi = s->mfHim[k];
            const float *xr = Xr + (size_t)k * MF_LANES, *xi = Xi + (size_t)k * MF_LANES;
            float *yr = Yr + (size_t)rev[k] * MF_LANES, *yi = Yi + (size_t)rev[k] * MF_LANES;
            for (int j = 0; j < MF_LANES; j++) {
                yr[j] = xr[j] * hr - xi[j] * hi;
                yi[j] = -(xr[j] * hi + xi[j] * hr);
            }
        }
        fft_rows(&s->mfPlan, Yr, Yi, MF_LANES, MF_LANES);
        for (int j = 0; j < nb; j++) {                          /* valid outputs, conjugated  */
            float *outR = s->cpiRe + (size_t)slot * s->nCells + (size_t)(b0 + j) * nR;
            float *outI = s->cpiIm + (size_t)slot * s->nCells + (size_t)(b0 + j) * nR;
            const float *yr = Yr + (size_t)(L - 1) * MF_LANES + j, *yi = Yi + (size_t)(L - 1) * MF_LANES + j;
            for (int r = 0; r < nR; r++) {
                outR[r] = yr[(size_t)r * MF_LANES];
                outI[r] = -yi[(size_t)r * MF_LANES];
            }
        }
    }
}

/* Valid convolution of every beam with the matched filter, written into the
 * CPI buffer slot: y[r] = sum_j h[j] * x[r + nTaps - 1 - j]. */
static void matched_filter(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot)
{
    if (s->mfFft) {
        matched_filter_fft(cfg, pc, s, slot);
        return;
    }
    const int n = pc->ddcCount, L = pc->nTaps, nR = s->nR;
    #pragma omp parallel for schedule(static)
    for (int b = 0; b < cfg->nBeams; b++) {
        const float *restrict xr = s->beamRe + (size_t)b * n, *restrict xi = s->beamIm + (size_t)b * n;
        float *restrict yr = s->cpiRe + (size_t)slot * s->nCells + (size_t)b * nR;
        float *restrict yi = s->cpiIm + (size_t)slot * s->nCells + (size_t)b * nR;
        memset(yr, 0, sizeof(float) * nR);
        memset(yi, 0, sizeof(float) * nR);
        for (int j = 0; j < L; j++) {
            const float hr = pc->tapsRe[j], hi = pc->tapsIm[j];
            const float *restrict ar = xr + (L - 1 - j), *restrict ai = xi + (L - 1 - j);
            for (int r = 0; r < nR; r++) {
                yr[r] += hr * ar[r] - hi * ai[r];
                yi[r] += hr * ai[r] + hi * ar[r];
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Slow time: MTI and Doppler                                               */
/* ------------------------------------------------------------------------ */

/* MTI output of the newest L pulses: sum_l c[stagger of the oldest][l] * x[oldest + l]. */
static void mti_newest(const RadarConfig *cfg, PulseState *s)
{
    const int L = cfg->nMtiTaps, N = s->nSlots;
    const int first = (s->last - (L - 1) + N) % N;               /* oldest pulse of the output */
    const double *c = cfg->mtiCoeffs + (size_t)(s->staggerIdx[first] - 1) * L;
    s->mtiLast = (s->mtiLast + 1) % s->nMtiOut;
    float *restrict yr = s->mtiRe + (size_t)s->mtiLast * s->nCells;
    float *restrict yi = s->mtiIm + (size_t)s->mtiLast * s->nCells;
    memset(yr, 0, sizeof(float) * s->nCells);
    memset(yi, 0, sizeof(float) * s->nCells);
    for (int l = 0; l < L; l++) {
        const float w = (float)c[l];
        const int slot = (first + l) % N;
        const float *restrict xr = s->cpiRe + (size_t)slot * s->nCells;
        const float *restrict xi = s->cpiIm + (size_t)slot * s->nCells;
        for (int k = 0; k < s->nCells; k++) {
            yr[k] += w * xr[k];
            yi[k] += w * xi[k];
        }
    }
    if (s->mtiCount < s->nMtiOut) s->mtiCount++;
}

/* Range-Doppler power map P[bin][cell] from the newest nRows*rowLen MTI outputs:
 * rows (fixed PRF: one row of consecutive outputs; staggered: de-interleaved by
 * stagger position), Hamming window, FFT, |X|^2 averaged over the rows, fftshift
 * and division by the noise gain of every bin. */
static void doppler(const Processor *p, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const FftPlan *plan = &p->fft;
    const int nfft = cfg->nfft, nRows = cfg->nRows, rowLen = cfg->rowLen;
    const int nUsed = nRows * rowLen;
    const int nBlocks = (s->nCells + CELL_BLOCK - 1) / CELL_BLOCK;
    const float scale = 1.0f / nRows;

    #pragma omp parallel num_threads(p->nThreads)
    {
        float *xr = p->scratch[thread_num()].xr;
        float *xi = p->scratch[thread_num()].xi;
        #pragma omp for schedule(static)
        for (int blk = 0; blk < nBlocks; blk++) {
            const int c0 = blk * CELL_BLOCK;
            const int nc = (s->nCells - c0 < CELL_BLOCK) ? s->nCells - c0 : CELL_BLOCK;
            for (int d = 0; d < nfft; d++)
                memset(s->P + (size_t)d * s->nCells + c0, 0, sizeof(float) * nc);
            for (int row = 0; row < nRows; row++) {
                /* load the row in bit-reversed order, windowed, zero padded */
                for (int j = 0; j < nfft; j++) {
                    float *dr = xr + (size_t)plan->bitrev[j] * CELL_BLOCK;
                    float *di = xi + (size_t)plan->bitrev[j] * CELL_BLOCK;
                    if (j >= rowLen) {
                        memset(dr, 0, sizeof(float) * nc);
                        memset(di, 0, sizeof(float) * nc);
                        continue;
                    }
                    const int q = j * nRows + row;                       /* chronological index */
                    const int slot = (s->mtiLast - (nUsed - 1) + q + 2 * s->nMtiOut) % s->nMtiOut;
                    const float w = cfg->dopWindow[j];
                    const float *sr = s->mtiRe + (size_t)slot * s->nCells + c0;
                    const float *si = s->mtiIm + (size_t)slot * s->nCells + c0;
                    for (int c = 0; c < nc; c++) {
                        dr[c] = w * sr[c];
                        di[c] = w * si[c];
                    }
                }
                fft_rows(plan, xr, xi, CELL_BLOCK, nc);
                for (int d = 0; d < nfft; d++) {                       /* fftshift + power */
                    const int sh = (d + nfft / 2) % nfft;
                    float *restrict pd = s->P + (size_t)sh * s->nCells + c0;
                    const float *restrict ar = xr + (size_t)d * CELL_BLOCK, *restrict ai = xi + (size_t)d * CELL_BLOCK;
                    for (int c = 0; c < nc; c++) pd[c] += ar[c] * ar[c] + ai[c] * ai[c];
                }
            }
            for (int d = 0; d < nfft; d++) {
                const float g = scale / cfg->noiseGain[d];
                float *restrict pd = s->P + (size_t)d * s->nCells + c0;
                for (int c = 0; c < nc; c++) pd[c] *= g;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* 2-D CA-CFAR                                                              */
/* ------------------------------------------------------------------------ */

/* y[r] = S[r+h+1] - S[r-h] with the indices clamped to 0 ... nR (zero padding);
 * S is the prefix sum of a range row (S[0] = 0). */
static void box_from_prefix(const double *S, float *y, int nR, int h)
{
    int r = 0;
    for (; r < nR && r < h; r++) y[r] = (float)(S[r + h + 1 < nR ? r + h + 1 : nR] - S[0]);
    for (; r + h + 1 <= nR; r++) y[r] = (float)(S[r + h + 1] - S[r - h]);
    for (; r < nR; r++) y[r] = (float)(S[nR] - S[r - h > 0 ? r - h : 0]);
}

/* Range box sums (outer and guard half widths) of four Doppler rows: prefix
 * sums in double, four independent chains so the additions overlap. */
static void range_box_sums4(const float *const x[4], double *S, float *const yO[4], float *const yG[4],
                            int nR, int hO, int hG)
{
    const int ld = nR + 1;
    double *S0 = S, *S1 = S + ld, *S2 = S + 2 * ld, *S3 = S + 3 * ld;
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    S0[0] = S1[0] = S2[0] = S3[0] = 0;
    for (int r = 0; r < nR; r++) {
        s0 += x[0][r];  s1 += x[1][r];  s2 += x[2][r];  s3 += x[3][r];
        S0[r + 1] = s0; S1[r + 1] = s1; S2[r + 1] = s2; S3[r + 1] = s3;
    }
    for (int i = 0; i < 4; i++) {
        box_from_prefix(S + (size_t)i * ld, yO[i], nR, hO);
        box_from_prefix(S + (size_t)i * ld, yG[i], nR, hG);
    }
}

/* 2-D CA-CFAR (as cfar_2d.m): noise = (sum over the outer window - sum over
 * the guard window) / nTrain; range zero padded, Doppler circular. One beam
 * at a time: range box sums of every Doppler bin, then running sums along
 * Doppler (double, no drift), fused with the threshold test. */
static void cfar(const Processor *p, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const int gR = cfg->guardRange, oR = gR + cfg->trainRange;
    const int gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
    const double alpha = cfg->alpha / cfg->nTrain;
    const int zeroBin = nfft / 2;
    const size_t nCells = (size_t)s->nCells;

    #pragma omp parallel num_threads(p->nThreads)
    {
        Scratch *w = &p->scratch[thread_num()];
        #pragma omp for schedule(dynamic)
        for (int b = 0; b < nB; b++) {
            for (int d0 = 0; d0 < nfft; d0 += 4) {
                const float *x[4];
                float *yO[4], *yG[4];
                for (int i = 0; i < 4; i++) {          /* rows past the end repeat the last one */
                    const int d = d0 + i < nfft ? d0 + i : nfft - 1;
                    x[i] = s->P + (size_t)d * nCells + (size_t)b * nR;
                    yO[i] = w->sumO + (size_t)d * nR;
                    yG[i] = w->sumG + (size_t)d * nR;
                }
                range_box_sums4(x, w->prefix, yO, yG, nR, oR, gR);
            }
            double *restrict aO = w->accO, *restrict aG = w->accG;
            for (int r = 0; r < nR; r++) aO[r] = aG[r] = 0.0;
            for (int k = -oD; k <= oD; k++) {                  /* window of bin 0 */
                const float *restrict o = w->sumO + (size_t)(((k % nfft) + nfft) % nfft) * nR;
                for (int r = 0; r < nR; r++) aO[r] += o[r];
            }
            for (int k = -gD; k <= gD; k++) {
                const float *restrict g = w->sumG + (size_t)(((k % nfft) + nfft) % nfft) * nR;
                for (int r = 0; r < nR; r++) aG[r] += g[r];
            }
            for (int d = 0; d < nfft; d++) {
                const float *restrict P = s->P + (size_t)d * nCells + (size_t)b * nR;
                unsigned char *restrict det = s->det + (size_t)d * nCells + (size_t)b * nR;
                if (abs(d - zeroBin) <= cfg->notchBins) {
                    memset(det, 0, nR);
                } else {
                    for (int r = 0; r < nR; r++) det[r] = P[r] > alpha * (aO[r] - aG[r]);
                    for (int r = 0; r < oR && r < nR; r++) det[r] = 0;          /* range edges */
                    for (int r = nR - oR > 0 ? nR - oR : 0; r < nR; r++) det[r] = 0;
                }
                /* slide the Doppler windows to bin d + 1 */
                const float *restrict addO = w->sumO + (size_t)((d + 1 + oD) % nfft) * nR;
                const float *restrict subO = w->sumO + (size_t)(((d - oD) % nfft + nfft) % nfft) * nR;
                const float *restrict addG = w->sumG + (size_t)((d + 1 + gD) % nfft) * nR;
                const float *restrict subG = w->sumG + (size_t)(((d - gD) % nfft + nfft) % nfft) * nR;
                for (int r = 0; r < nR; r++) {
                    aO[r] += (double)addO[r] - subO[r];
                    aG[r] += (double)addG[r] - subG[r];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Pruning and detection reports                                            */
/* ------------------------------------------------------------------------ */

typedef struct { int b, r, d; float pw; } Cell;

/* MATLAB's find() order: Doppler, then range, then beam. */
static int cell_order(const void *a, const void *b)
{
    const Cell *x = a, *y = b;
    if (x->d != y->d) return x->d - y->d;
    if (x->r != y->r) return x->r - y->r;
    return x->b - y->b;
}

static float P_at(const PulseState *s, int d, int b, int r)
{
    return s->P[(size_t)d * s->nCells + (size_t)b * s->nR + r];
}

/* k-th smallest value (0-based) of v[0..n-1]; v is reordered. */
static float select_kth(float *v, int n, int k)
{
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        float pivot = v[(lo + hi) / 2];
        int i = lo, j = hi;
        while (i <= j) {
            while (v[i] < pivot) i++;
            while (v[j] > pivot) j--;
            if (i <= j) { float t = v[i]; v[i] = v[j]; v[j] = t; i++; j--; }
        }
        if (k <= j) hi = j; else if (k >= i) lo = i; else return v[k];
    }
    return v[k];
}

/* Noise floor of a beam: median of its range-Doppler map / (median/mean of the noise). */
static double beam_noise_floor(const RadarConfig *cfg, const PulseState *s, int b, float *work)
{
    const int n = cfg->nfft * s->nR;
    for (int d = 0; d < cfg->nfft; d++)
        memcpy(work + (size_t)d * s->nR, s->P + (size_t)d * s->nCells + (size_t)b * s->nR, sizeof(float) * s->nR);
    double med;
    if (n % 2) {
        med = select_kth(work, n, n / 2);
    } else {                                                    /* mean of the two middle values */
        float hi = select_kth(work, n, n / 2);
        float lo = work[0];
        for (int i = 1; i < n / 2; i++) if (work[i] > lo) lo = work[i];   /* max of the lower half */
        med = 0.5 * ((double)lo + hi);
    }
    return med / cfg->medianToMean;
}

static void add_report(PulseState *s, const Detection *d)
{
    if (s->nReports == s->capReports) {
        s->capReports = s->capReports ? 2 * s->capReports : 64;
        s->reports = xrealloc(s->reports, sizeof(Detection) * s->capReports);
    }
    s->reports[s->nReports++] = *d;
}

/* Keeps the detections that are local peaks (3x3x3 beam x range x Doppler) and
 * not sidelobes of a much stronger peak, and converts them to reports. */
static void prune_and_report(const RadarConfig *cfg, int pulse, PulseState *s)
{
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const PulseConfig *pc = &cfg->pulse[pulse];
    Cell *cells = NULL;
    int n = 0, cap = 0;

    /* 1) local peaks (scanned in memory order, then sorted into the order of
     *    MATLAB's find() on [beam x range x Doppler]) */
    for (int d = 0; d < nfft; d++) {
        for (int b = 0; b < nB; b++) {
            const unsigned char *detRow = s->det + (size_t)d * s->nCells + (size_t)b * nR;
            const unsigned char *hit = memchr(detRow, 1, nR);
            for (int r = hit ? (int)(hit - detRow) : nR; r < nR; r++) {
                if (!detRow[r]) continue;
                const float v = P_at(s, d, b, r);
                int peak = 1;
                for (int dd = -1; dd <= 1 && peak; dd++) {
                    const int d2 = (d + dd + nfft) % nfft;
                    for (int bb = b - 1; bb <= b + 1 && peak; bb++) {
                        if (bb < 0 || bb >= nB) continue;
                        for (int rr = r - 1; rr <= r + 1; rr++) {
                            if (rr < 0 || rr >= nR) continue;
                            if (P_at(s, d2, bb, rr) > v) { peak = 0; break; }
                        }
                    }
                }
                if (!peak) continue;
                if (n == cap) {
                    cap = cap ? 2 * cap : 64;
                    cells = xrealloc(cells, sizeof(Cell) * cap);
                }
                cells[n++] = (Cell){ b, r, d, v };
            }
        }
    }
    if (n > 1) qsort(cells, n, sizeof(Cell), cell_order);

    /* 2) sidelobe blanking (range, elevation and Doppler sidelobes) */
    const double ratio = pow(10.0, cfg->sidelobeMargin_dB / 10.0);
    const int pulseCells = pc->nTaps;
    unsigned char *keep = xcalloc(n, 1);
    for (int i = 0; i < n; i++) {
        keep[i] = 1;
        for (int j = 0; j < n && keep[i]; j++) {
            if (!((double)cells[j].pw > ratio * cells[i].pw)) continue;
            int dDop = abs(cells[i].d - cells[j].d);
            if (nfft - dDop < dDop) dDop = nfft - dDop;
            const int sameD = dDop <= 1, sameB = abs(cells[i].b - cells[j].b) <= 1;
            const int dR = abs(cells[i].r - cells[j].r), sameR = dR <= 1;
            if ((sameB && sameD && dR <= pulseCells) || (sameR && sameD) || (sameR && sameB)) keep[i] = 0;
        }
    }

    /* 3) reports: range, velocity, elevation, CPI centre azimuth/time, SNR */
    s->nReports = 0;
    double sinSum = 0, cosSum = 0, tSum = 0;
    for (int k = 0; k < s->nSlots; k++) {
        sinSum += sin(s->az[k] * PI / 180);
        cosSum += cos(s->az[k] * PI / 180);
        tSum += s->time[k];
    }
    double azCenter = atan2(sinSum / s->nSlots, cosSum / s->nSlots) * 180 / PI;
    azCenter = fmod(azCenter + 360.0, 360.0);
    double *floors = NULL;
    float *work = NULL;
    for (int i = 0; i < n; i++) {
        if (!keep[i]) continue;
        const Cell *c = &cells[i];
        if (!floors) {
            floors = xcalloc(nB, sizeof(double));
            for (int b = 0; b < nB; b++) floors[b] = -1;
            work = xcalloc((size_t)nfft * nR, sizeof(float));
        }
        if (floors[c->b] < 0) floors[c->b] = beam_noise_floor(cfg, s, c->b, work);
        Detection d;
        memset(&d, 0, sizeof d);
        d.time_s = tSum / s->nSlots;
        d.az_deg = azCenter;
        d.range_m = pc->range_m[c->r];
        d.vr_mps = -cfg->dopFreq_Hz[c->d] * pc->lambda_m / 2;
        d.el_deg = NAN;
        if (cfg->dbfEnabled) {
            double u = cfg->beamU[c->b];
            if (c->b > 0 && c->b < nB - 1) {                   /* parabola through the beam powers [dB] */
                const double p0 = 10 * log10(P_at(s, c->d, c->b - 1, c->r));
                const double p1 = 10 * log10(P_at(s, c->d, c->b, c->r));
                const double p2 = 10 * log10(P_at(s, c->d, c->b + 1, c->r));
                const double den = p0 - 2 * p1 + p2;
                const double pos = c->b + 0.5 * (p0 - p2) / den;       /* fractional beam index */
                if (den < 0 && isfinite(pos)) {                        /* (a zero power gives -Inf dB) */
                    const int i0 = (int)floor(pos);
                    const double f = pos - i0;
                    u = (i0 + 1 < nB) ? cfg->beamU[i0] * (1 - f) + cfg->beamU[i0 + 1] * f : cfg->beamU[i0];
                }
            }
            d.el_deg = cfg->tilt_deg + asin(u) * 180 / PI;
        }
        double ratioSnr = c->pw / floors[c->b] - 1;
        d.snr_dB = 10 * log10(ratioSnr > DBL_EPSILON ? ratioSnr : DBL_EPSILON);
        d.power = c->pw;
        d.beam = c->b;  d.rIdx = c->r;  d.dIdx = c->d;
        d.pulse = pulse;
        d.hits = 0;
        d.pulseNumber = s->newestPulse;
        add_report(s, &d);
    }
    free(floors);
    free(work);
    free(keep);
    free(cells);
}

/* ------------------------------------------------------------------------ */
/* One PRI                                                                  */
/* ------------------------------------------------------------------------ */

static double elapsed_ms(double *t0)
{
    double t = now_seconds(), dt = (t - *t0) * 1e3;
    *t0 = t;
    return dt;
}

static void plot_callback(void *user, const Detection *plot)
{
    const ProcessorCallbacks *cb = user;
    if (cb && cb->plot) cb->plot(cb->user, plot);
}

/* The PRI has the layout of the configuration: channels, DDC samples, the DDC
 * samples of both sub-pulses (header, 1-based) and the DDC rate. */
static int pri_matches_config(const RadarConfig *cfg, const PriData *pri)
{
    if (pri->nCh != cfg->nCh || pri->nSamples != cfg->nDdc || !pri->iq) return 0;
    const uint32_t *idx[N_PULSES] = { pri->shortIdx, pri->longIdx };
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        if (idx[i][0] != (uint32_t)pc->ddcFirst + 1 || idx[i][1] != (uint32_t)(pc->ddcFirst + pc->ddcCount))
            return 0;
    }
    return fabs(pri->ddcFs_Hz - cfg->fsOut_Hz) <= 1e-3 * cfg->fsOut_Hz;
}

void processor_flush(Processor *p, const ProcessorCallbacks *cb)
{
    for (int i = 0; i < N_PULSES; i++)
        extractor_update(&p->ps[i].ext, NULL, 0, 1, plot_callback, (void *)cb);
}

void processor_process(Processor *p, const PriData *pri, const ProcessorCallbacks *cb)
{
    const RadarConfig *cfg = p->cfg;
    if (!pri_matches_config(cfg, pri)) {                    /* made for another configuration */
        p->timing.rejected++;
        return;
    }

    /* a gap in the pulse numbers (skipped sector): restart the CPIs, close the plots */
    if (p->havePulse && pri->pulseNumber != p->lastPulse + 1) {
        processor_flush(p, cb);
        for (int i = 0; i < N_PULSES; i++) reset_buffers(&p->ps[i]);
    }
    p->lastPulse = pri->pulseNumber;
    p->havePulse = 1;
    if (cb && cb->antenna) cb->antenna(cb->user, pri->az_deg, pri->time_s, pri->scanNumber);

    const double tStart = now_seconds();
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        double t0 = now_seconds();

        /* fast time into the next CPI slot */
        load_segment(cfg, pc, pri, s);
        beamform(cfg, pc, s);
        p->timing.dbf_ms += elapsed_ms(&t0);
        const int slot = (s->last + 1) % s->nSlots;
        matched_filter(cfg, pc, s, slot);
        s->last = slot;
        s->staggerIdx[slot] = pri->staggerIndex >= 1 && pri->staggerIndex <= cfg->nStagger ? pri->staggerIndex : 1;
        s->az[slot] = pri->az_deg;
        s->time[slot] = pri->time_s;
        s->newestPulse = pri->pulseNumber;
        if (s->count < s->nSlots) s->count++;
        s->nPushed = (s->nPushed + 1) % cfg->updateEvery;  /* pulses since the last update (mod N) */
        p->timing.mf_ms += elapsed_ms(&t0);

        /* slow time */
        if (s->count >= cfg->nMtiTaps) mti_newest(cfg, s);
        p->timing.mti_ms += elapsed_ms(&t0);
        if (s->count < s->nSlots || s->nPushed != 0) continue;
        doppler(p, s);
        p->timing.doppler_ms += elapsed_ms(&t0);
        cfar(p, s);
        prune_and_report(cfg, i, s);
        p->timing.cfar_ms += elapsed_ms(&t0);
        if (cb && cb->detections) cb->detections(cb->user, i, s->reports, s->nReports);
        extractor_update(&s->ext, s->reports, s->nReports, 0, plot_callback, (void *)cb);
        p->timing.extract_ms += elapsed_ms(&t0);
        if (i == 0) p->timing.updates++;
    }
    p->timing.total_ms += (now_seconds() - tStart) * 1e3;
    p->timing.pris++;
}
