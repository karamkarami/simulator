/*
 * fft.c - radix-2 decimation-in-time FFT of many signals at once.
 */
#include "fft.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

int fft_plan_init(FftPlan *plan, int n)
{
    plan->n = n;
    plan->bitrev = malloc(sizeof(int) * n);
    plan->twRe = malloc(sizeof(float) * (n / 2 + 1));
    plan->twIm = malloc(sizeof(float) * (n / 2 + 1));
    if (!plan->bitrev || !plan->twRe || !plan->twIm || n < 1 || (n & (n - 1))) return -1;
    int bits = 0;
    while ((1 << bits) < n) bits++;
    for (int m = 0; m < n; m++) {
        int r = 0;
        for (int b = 0; b < bits; b++)
            if (m & (1 << b)) r |= 1 << (bits - 1 - b);
        plan->bitrev[m] = r;
    }
    for (int k = 0; k < n / 2; k++) {
        double a = -2.0 * PI * k / n;
        plan->twRe[k] = (float)cos(a);
        plan->twIm[k] = (float)sin(a);
    }
    return 0;
}

void fft_plan_free(FftPlan *plan)
{
    free(plan->bitrev);
    free(plan->twRe);
    free(plan->twIm);
    plan->bitrev = NULL;
    plan->twRe = plan->twIm = NULL;
}

/* The butterflies; with a constant nSignals the compiler unrolls the inner loop. */
#if defined(__GNUC__)
__attribute__((always_inline))
#endif
static inline void fft_core(const FftPlan *plan, float *re, float *im, size_t stride, int nSignals)
{
    const int n = plan->n;
    for (int size = 2; size <= n; size *= 2) {
        const int half = size / 2;
        const int step = n / size;                     /* twiddle index step */
        for (int start = 0; start < n; start += size) {
            for (int k = 0; k < half; k++) {
                const float wr = plan->twRe[k * step];
                const float wi = plan->twIm[k * step];
                float *restrict ar = re + (size_t)(start + k) * stride;
                float *restrict ai = im + (size_t)(start + k) * stride;
                float *restrict br = re + (size_t)(start + k + half) * stride;
                float *restrict bi = im + (size_t)(start + k + half) * stride;
                for (int c = 0; c < nSignals; c++) {   /* one butterfly per signal */
                    const float tr = wr * br[c] - wi * bi[c];
                    const float ti = wr * bi[c] + wi * br[c];
                    br[c] = ar[c] - tr;
                    bi[c] = ai[c] - ti;
                    ar[c] += tr;
                    ai[c] += ti;
                }
            }
        }
    }
}

void fft_rows(const FftPlan *plan, float *re, float *im, size_t stride, int nSignals)
{
    switch (nSignals) {                                /* common widths: unrolled */
    case 24: fft_core(plan, re, im, stride, 24); break;
    case 32: fft_core(plan, re, im, stride, 32); break;
    case 64: fft_core(plan, re, im, stride, 64); break;
    default: fft_core(plan, re, im, stride, nSignals);
    }
}
