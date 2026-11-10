/* C188 (Claude Opus 5.5): matrix-free momentum-sector Lanczos for a translation-invariant qubit ring
   H = sum_i T^i h T^-i, h a 4x4 complex bond on sites (i, i+1 mod L).
   Conventions copied from C155's RingTable (est_c155.py): site i is bit (L-1-i); local index of a bond
   state is 2*b_i + b_{i+1}; <y|H|x> = h[out][in]; T: cur -> (cur>>1)|((cur&1)<<(L-1)); x = T^shift rep;
   sector-m matrix element M[ry][R] = h[out][in] e^{ik shift(y)} sqrt(per_R/per_ry).
   We apply the Hermitian M row-wise: (Mv)[R'] = sum conj(M[ry][R']) v[ry] over the images y of R'.
   Lanczos without reorthogonalization; alpha/beta are written out and the spectrum (with Cullum-Willoughby
   ghost removal) is computed by ringk_post.py.
   Usage: ringk L m niter seed hfile out
     hfile: 32 numbers, Re and Im of h[out][in] row-major (out*4+in), as "re im" pairs.
   Output: text, first line "L m n", then niter lines "alpha beta". */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef double complex cplx;

static inline uint64_t rot(uint64_t c, int L, uint64_t mask) { return ((c >> 1) | ((c & 1ULL) << (L - 1))) & mask; }

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: ringk L m niter seed hfile out\n"); return 1; }
    int L = atoi(argv[1]); int m = atoi(argv[2]); int niter = atoi(argv[3]); unsigned seed = (unsigned)atoi(argv[4]);
    cplx h[16];
    FILE *fh = fopen(argv[5], "r");
    for (int i = 0; i < 16; i++) { double re, im; if (fscanf(fh, "%lf %lf", &re, &im) != 2) { fprintf(stderr, "bad h\n"); return 1; } h[i] = re + I * im; }
    fclose(fh);
    uint64_t N = 1ULL << L, mask = N - 1;
    double k = 2.0 * M_PI * m / L;
    /* pass A: flag reps of this sector, store period */
    uint8_t *per = (uint8_t *)malloc(N);
    #pragma omp parallel for schedule(static)
    for (int64_t x = 0; x < (int64_t)N; x++) {
        uint64_t c = (uint64_t)x; int isrep = 1, p = L;
        for (int j = 1; j < L; j++) {
            c = rot(c, L, mask);
            if (c < (uint64_t)x) { isrep = 0; break; }
            if (c == (uint64_t)x) { p = j; break; }
        }
        per[x] = (isrep && ((m * p) % L == 0)) ? (uint8_t)p : 0;
    }
    /* prefix: positions */
    int64_t n = 0;
    for (uint64_t x = 0; x < N; x++) if (per[x]) n++;
    uint32_t *reps = (uint32_t *)malloc(n * sizeof(uint32_t));
    uint8_t *pers = (uint8_t *)malloc(n);
    int32_t *tab = (int32_t *)malloc(N * sizeof(int32_t));
    { int64_t q = 0; for (uint64_t x = 0; x < N; x++) { if (per[x]) { reps[q] = (uint32_t)x; pers[q] = per[x]; tab[x] = (int32_t)(q * 32); q++; } else tab[x] = -2; } }
    free(per);
    /* pass B: every x -> pos(rep)*32 + shift, or -1 */
    #pragma omp parallel for schedule(static)
    for (int64_t x = 0; x < (int64_t)N; x++) {
        if (tab[x] >= 0) continue;   /* a rep of this sector, shift 0 */
        uint64_t c = (uint64_t)x, best = c; int jmin = 0;
        for (int j = 1; j < L; j++) { c = rot(c, L, mask); if (c < best) { best = c; jmin = j; } }
        int32_t tb = (best == (uint64_t)x) ? -2 : tab[best];
        if (tb == -2 || tb == -1) { tab[x] = -1; continue; }
        /* tb was set in the serial loop (rep of the sector, value q*32) and is never rewritten */
        int shift = (L - jmin) % L;
        tab[x] = (tb & ~31) + shift;
    }
    #pragma omp parallel for schedule(static)
    for (int64_t x = 0; x < (int64_t)N; x++) if (tab[x] == -2) tab[x] = -1;
    fprintf(stderr, "L=%d m=%d sector dim %lld\n", L, m, (long long)n);
    cplx ph[64]; for (int s = 0; s < L; s++) ph[s] = cexp(I * k * s);
    double *sq = (double *)malloc(64 * 64 * sizeof(double));
    for (int a = 1; a <= L; a++) for (int b = 1; b <= L; b++) sq[a * 64 + b] = sqrt((double)a / (double)b);
    /* nonzero entries of h by input column */
    int nout[4]; int outs[4][4]; cplx hv[4][4];
    for (int in = 0; in < 4; in++) { nout[in] = 0; for (int o = 0; o < 4; o++) if (cabs(h[o * 4 + in]) > 1e-14) { outs[in][nout[in]] = o; hv[in][nout[in]] = conj(h[o * 4 + in]); nout[in]++; } }
    cplx *v = (cplx *)malloc(n * sizeof(cplx)), *vp = (cplx *)malloc(n * sizeof(cplx)), *w = (cplx *)malloc(n * sizeof(cplx));
    srand(seed);
    double nrm = 0;
    for (int64_t i = 0; i < n; i++) { double a = rand() / (double)RAND_MAX - 0.5, b = rand() / (double)RAND_MAX - 0.5; v[i] = a + I * b; nrm += a * a + b * b; }
    nrm = sqrt(nrm);
    for (int64_t i = 0; i < n; i++) { v[i] /= nrm; vp[i] = 0; }
    FILE *fo = fopen(argv[6], "w");
    fprintf(fo, "%d %d %lld\n", L, m, (long long)n);
    double beta_prev = 0;
    for (int it = 0; it < niter; it++) {
        /* w = M v */
        #pragma omp parallel for schedule(dynamic, 4096)
        for (int64_t r = 0; r < n; r++) {
            uint64_t R = reps[r]; int pR = pers[r]; cplx acc = 0;
            for (int i = 0; i < L; i++) {
                int b0 = L - 1 - i, b1 = L - 1 - ((i + 1) % L);
                int in = (int)(((R >> b0) & 1) << 1 | ((R >> b1) & 1));
                uint64_t base = R & ~(1ULL << b0) & ~(1ULL << b1);
                for (int t = 0; t < nout[in]; t++) {
                    int o = outs[in][t];
                    uint64_t y = base | ((uint64_t)((o >> 1) & 1) << b0) | ((uint64_t)(o & 1) << b1);
                    int32_t tb = tab[y]; if (tb < 0) continue;
                    int64_t q = tb >> 5; int s = tb & 31;
                    acc += hv[in][t] * conj(ph[s]) * sq[pR * 64 + pers[q]] * v[q];
                }
            }
            w[r] = acc;
        }
        double alpha = 0;
        #pragma omp parallel for reduction(+:alpha)
        for (int64_t i = 0; i < n; i++) alpha += creal(conj(v[i]) * w[i]);
        double bn = 0;
        #pragma omp parallel for reduction(+:bn)
        for (int64_t i = 0; i < n; i++) { w[i] -= alpha * v[i] + beta_prev * vp[i]; bn += creal(conj(w[i]) * w[i]); }
        bn = sqrt(bn);
        fprintf(fo, "%.17g %.17g\n", alpha, bn); fflush(fo);
        if (bn < 1e-13) break;
        #pragma omp parallel for
        for (int64_t i = 0; i < n; i++) { vp[i] = v[i]; v[i] = w[i] / bn; }
        beta_prev = bn;
    }
    fclose(fo);
    return 0;
}
