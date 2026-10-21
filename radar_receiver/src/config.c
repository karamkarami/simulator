/*
 * config.c - loads the receiver configuration written by export_receiver_config.m.
 *
 * File: 8 bytes "RADRXCFG", uint32 version, uint32 number of records; every
 * record: char name[32], uint32 type (0 real, 1 complex), uint32 rows,
 * uint32 cols, rows*cols doubles (complex: re, im pairs), column-major.
 */
#include "radar_rx.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char     name[33];
    uint32_t isComplex, rows, cols;
    double  *data;          /* rows*cols values (complex: 2*rows*cols) */
} Record;

typedef struct {
    Record *rec;
    int     n;
    char   *err;
    size_t  errLen;
    int     failed;
} RecordSet;

static void set_error(RecordSet *rs, const char *fmt, ...)
{
    if (rs->failed) return;                       /* keep the first error */
    rs->failed = 1;
    va_list ap;
    va_start(ap, fmt);
    if (rs->err && rs->errLen) vsnprintf(rs->err, rs->errLen, fmt, ap);
    va_end(ap);
}

static const Record *find(RecordSet *rs, const char *name)
{
    for (int i = 0; i < rs->n; i++)
        if (strcmp(rs->rec[i].name, name) == 0) return &rs->rec[i];
    set_error(rs, "configuration value '%s' is missing", name);
    return NULL;
}

static double get_scalar(RecordSet *rs, const char *name)
{
    const Record *r = find(rs, name);
    return (r && r->rows * r->cols >= 1) ? r->data[0] : 0.0;
}

static int get_int(RecordSet *rs, const char *name)
{
    return (int)lround(get_scalar(rs, name));
}

/* Real vector of n values (n = 0: any length, returned in *nOut). */
static double *get_doubles(RecordSet *rs, const char *name, int n, int *nOut)
{
    const Record *r = find(rs, name);
    if (!r) return NULL;
    int count = (int)(r->rows * r->cols);
    if (r->isComplex || (n > 0 && count != n)) {
        set_error(rs, "configuration value '%s' has %d values, expected %d", name, count, n);
        return NULL;
    }
    double *v = malloc(sizeof(double) * (count > 0 ? count : 1));
    if (!v) { set_error(rs, "out of memory"); return NULL; }
    memcpy(v, r->data, sizeof(double) * count);
    if (nOut) *nOut = count;
    return v;
}

static float *get_floats(RecordSet *rs, const char *name, int n)
{
    double *d = get_doubles(rs, name, n, NULL);
    if (!d) return NULL;
    float *f = malloc(sizeof(float) * (n > 0 ? n : 1));
    if (!f) { free(d); set_error(rs, "out of memory"); return NULL; }
    for (int i = 0; i < n; i++) f[i] = (float)d[i];
    free(d);
    return f;
}

/* Complex values as separate real / imaginary float arrays. rowMajor: convert
 * a [rows x cols] MATLAB matrix (column-major) to row-major order. */
static int get_complex(RecordSet *rs, const char *name, int rows, int cols, int rowMajor,
                       float **re, float **im)
{
    const Record *r = find(rs, name);
    if (!r) return -1;
    if (!r->isComplex || (int)r->rows != rows || (int)r->cols != cols) {
        set_error(rs, "configuration value '%s' should be complex %d x %d", name, rows, cols);
        return -1;
    }
    *re = malloc(sizeof(float) * rows * cols);
    *im = malloc(sizeof(float) * rows * cols);
    if (!*re || !*im) { set_error(rs, "out of memory"); return -1; }
    for (int c = 0; c < cols; c++) {
        for (int rr = 0; rr < rows; rr++) {
            int src = c * rows + rr;                       /* column-major */
            int dst = rowMajor ? rr * cols + c : src;
            (*re)[dst] = (float)r->data[2 * src];
            (*im)[dst] = (float)r->data[2 * src + 1];
        }
    }
    return 0;
}

static int read_records(FILE *f, RecordSet *rs)
{
    char magic[8];
    uint32_t head[2];
    if (fread(magic, 1, 8, f) != 8 || memcmp(magic, "RADRXCFG", 8) != 0 ||
        fread(head, sizeof(uint32_t), 2, f) != 2 || head[0] != 1) {
        set_error(rs, "not a receiver configuration file (version 1)");
        return -1;
    }
    rs->n = (int)head[1];
    rs->rec = calloc(rs->n > 0 ? rs->n : 1, sizeof(Record));
    if (!rs->rec) { set_error(rs, "out of memory"); return -1; }
    for (int i = 0; i < rs->n; i++) {
        Record *r = &rs->rec[i];
        uint32_t dims[3];
        if (fread(r->name, 1, 32, f) != 32 || fread(dims, sizeof(uint32_t), 3, f) != 3) {
            set_error(rs, "truncated configuration file");
            return -1;
        }
        r->name[32] = '\0';
        r->isComplex = dims[0];
        r->rows = dims[1];
        r->cols = dims[2];
        size_t count = (size_t)r->rows * r->cols * (r->isComplex ? 2 : 1);
        r->data = malloc(sizeof(double) * (count > 0 ? count : 1));
        if (!r->data || fread(r->data, sizeof(double), count, f) != count) {
            set_error(rs, "truncated configuration file (record '%s')", r->name);
            return -1;
        }
    }
    return 0;
}

static void free_records(RecordSet *rs)
{
    for (int i = 0; i < rs->n; i++) free(rs->rec[i].data);
    free(rs->rec);
}

int radar_config_load(const char *path, RadarConfig *cfg, char *err, size_t errLen)
{
    RecordSet rs = { NULL, 0, err, errLen, 0 };
    memset(cfg, 0, sizeof(*cfg));
    FILE *f = fopen(path, "rb");
    if (!f) {
        if (err && errLen) snprintf(err, errLen, "cannot open %s", path);
        return -1;
    }
    int ok = read_records(f, &rs) == 0;
    fclose(f);

    if (ok) {
        cfg->nCh             = get_int(&rs, "nCh");
        cfg->nBeams          = get_int(&rs, "nBeams");
        cfg->dbfEnabled      = get_int(&rs, "dbfEnabled");
        if (cfg->dbfEnabled)
            get_complex(&rs, "W", cfg->nBeams, cfg->nCh, 1, &cfg->wRe, &cfg->wIm);
        cfg->beamU           = get_doubles(&rs, "beamU", cfg->nBeams, NULL);
        cfg->tilt_deg        = get_scalar(&rs, "tilt_deg");
        cfg->height_m        = get_scalar(&rs, "height_m");
        cfg->rate_dps        = get_scalar(&rs, "rate_dps");
        cfg->azBeamwidth_deg = get_scalar(&rs, "azBeamwidth_deg");
        cfg->maxRange_m      = get_scalar(&rs, "maxRange_m");
        cfg->nDdc            = get_int(&rs, "nDdc");
        cfg->fsOut_Hz        = get_scalar(&rs, "fsOut_Hz");

        static const char *names[N_PULSES] = { "short", "long" };
        for (int p = 0; p < N_PULSES; p++) {
            PulseConfig *pc = &cfg->pulse[p];
            char key[48];
            snprintf(key, sizeof key, "%s_ddcFirst", names[p]);  pc->ddcFirst = get_int(&rs, key);
            snprintf(key, sizeof key, "%s_ddcCount", names[p]);  pc->ddcCount = get_int(&rs, key);
            snprintf(key, sizeof key, "%s_nTaps", names[p]);     pc->nTaps    = get_int(&rs, key);
            snprintf(key, sizeof key, "%s_nRange", names[p]);    pc->nRange   = get_int(&rs, key);
            snprintf(key, sizeof key, "%s_taps", names[p]);
            get_complex(&rs, key, 1, pc->nTaps, 0, &pc->tapsRe, &pc->tapsIm);
            snprintf(key, sizeof key, "%s_range_m", names[p]);  pc->range_m = get_doubles(&rs, key, pc->nRange, NULL);
            snprintf(key, sizeof key, "%s_delay_s", names[p]);  pc->delay_s = get_doubles(&rs, key, pc->nRange, NULL);
            snprintf(key, sizeof key, "%s_lambda_m", names[p]); pc->lambda_m = get_scalar(&rs, key);
            if (!rs.failed && (pc->ddcCount != pc->nRange + pc->nTaps - 1 ||
                               pc->ddcFirst < 0 || pc->ddcFirst + pc->ddcCount > cfg->nDdc))
                set_error(&rs, "inconsistent DDC segment of the %s pulse", names[p]);
        }

        cfg->cpiLength  = get_int(&rs, "cpiLength");
        cfg->nStagger   = get_int(&rs, "nStagger");
        cfg->staggered  = get_int(&rs, "staggered");
        {
            const Record *r = find(&rs, "mtiCoeffs");
            if (r && !r->isComplex && (int)r->rows == cfg->nStagger) {
                cfg->nMtiTaps = (int)r->cols;
                cfg->mtiCoeffs = malloc(sizeof(double) * r->rows * r->cols);
                for (uint32_t k = 0; k < r->rows; k++)          /* to row-major */
                    for (uint32_t l = 0; l < r->cols; l++)
                        cfg->mtiCoeffs[k * r->cols + l] = r->data[l * r->rows + k];
            } else if (r) {
                set_error(&rs, "mtiCoeffs must have nStagger rows");
            }
        }
        cfg->nRows        = get_int(&rs, "nRows");
        cfg->rowLen       = get_int(&rs, "rowLen");
        cfg->nfft         = get_int(&rs, "nfft");
        cfg->dopWindow    = get_floats(&rs, "dopWindow", cfg->rowLen);
        cfg->noiseGain    = get_floats(&rs, "noiseGain", cfg->nfft);
        cfg->dopFreq_Hz   = get_doubles(&rs, "dopFreq_Hz", cfg->nfft, NULL);
        cfg->notchBins    = get_int(&rs, "notchBins");
        cfg->medianToMean = get_scalar(&rs, "medianToMean");
        cfg->guardRange   = get_int(&rs, "guardRange");
        cfg->trainRange   = get_int(&rs, "trainRange");
        cfg->guardDoppler = get_int(&rs, "guardDoppler");
        cfg->trainDoppler = get_int(&rs, "trainDoppler");
        cfg->nTrain       = get_int(&rs, "nTrain");
        cfg->alpha        = get_scalar(&rs, "alpha");
        cfg->sidelobeMargin_dB = get_scalar(&rs, "sidelobeMargin_dB");
        cfg->minHits      = get_int(&rs, "minHits");
        cfg->maxMiss      = get_int(&rs, "maxMiss");
        cfg->updateEvery  = get_int(&rs, "updateEvery");
        cfg->minExtent_deg = get_scalar(&rs, "minExtent_deg");
        cfg->meanPri_s    = get_scalar(&rs, "meanPri_s");

        if (!rs.failed && (cfg->updateEvery < 1 || cfg->minHits < 1 || cfg->cpiLength < 1 || cfg->nfft < 1 ||
                           cfg->nRows < 1 || cfg->nBeams < 1 || cfg->nCh < 1))
            set_error(&rs, "invalid sizes in the configuration");
        if (!rs.failed && (cfg->nfft & (cfg->nfft - 1)) != 0)
            set_error(&rs, "the Doppler FFT length %d is not a power of 2", cfg->nfft);
        if (!rs.failed && cfg->nRows * cfg->rowLen + cfg->nMtiTaps - 1 > cfg->cpiLength)
            set_error(&rs, "the CPI buffer is shorter than the Doppler processing needs");
    }
    free_records(&rs);
    if (!ok || rs.failed) {
        radar_config_free(cfg);
        return -1;
    }
    return 0;
}

void radar_config_free(RadarConfig *cfg)
{
    free(cfg->wRe);  free(cfg->wIm);  free(cfg->beamU);
    for (int p = 0; p < N_PULSES; p++) {
        PulseConfig *pc = &cfg->pulse[p];
        free(pc->tapsRe);  free(pc->tapsIm);  free(pc->range_m);  free(pc->delay_s);
    }
    free(cfg->mtiCoeffs);  free(cfg->dopWindow);  free(cfg->noiseGain);  free(cfg->dopFreq_Hz);
    memset(cfg, 0, sizeof(*cfg));
}
