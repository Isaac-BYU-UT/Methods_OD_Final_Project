/*
 * MEX implementation of funcGottliebNorm.m
 * Jack Li, jack.li@utexas.edu
 * April 2026
 *
 * Calculates acceleration due to Earth's gravity using NASA's
 * normalised Gottlieb algorithm.
 *
 * Usage (identical to MATLAB version):
 *   aGravity = funcGottliebNormMEX(mu, re, xin, C, S, nax, mmax, rnp)
 *
 * Inputs:
 *   mu   - gravitational parameter (km^3/s^2), scalar
 *   re   - Earth radius (km), scalar
 *   xin  - position vector [3x1] in km (only first 3 elements used)
 *   C    - C coefficients matrix [(nax+1) x (nax+1)]
 *   S    - S coefficients matrix [(nax+1) x (nax+1)]
 *   nax  - maximum degree n, scalar
 *   mmax - maximum order m, scalar
 *   rnp  - rotation matrix [3x3] (ECEF to ECI, typically eye(3))
 *
 * Output:
 *   aGravity - gravity acceleration vector [3x1] in km/s^2
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Helper: column-major index for matrix stored as 1D array */
/* MATLAB stores matrices column-major: A(i,j) = A[j*nrows + i] */
#define IDX(i, j, nrows) ((j)*(nrows) + (i))

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ---- Input validation ---- */
    if (nrhs != 8)
        mexErrMsgIdAndTxt("funcGottliebNorm:nrhs",
                          "Eight inputs required: mu, re, xin, C, S, nax, mmax, rnp");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("funcGottliebNorm:nlhs",
                          "One output required.");

    /* ---- Read scalar inputs ---- */
    double mu   = mxGetScalar(prhs[0]);
    double re   = mxGetScalar(prhs[1]);
    int    nax  = (int)mxGetScalar(prhs[5]);
    int    mmax = (int)mxGetScalar(prhs[6]);

    /* ---- Read xin (position vector, use first 3 elements) ---- */
    double *xin_ptr = mxGetPr(prhs[2]);

    /* ---- Read C and S matrices ---- */
    double *C_mat = mxGetPr(prhs[3]);
    double *S_mat = mxGetPr(prhs[4]);
    int     CS_rows = (int)mxGetM(prhs[3]);   /* number of rows in C */

    /* ---- Read rnp rotation matrix (3x3) ---- */
    double *rnp = mxGetPr(prhs[7]);

    /* ---- Allocate working arrays ---- */
    /* All indexed 1..nax+1 in MATLAB — we allocate nax+2 for safety */
    int sz = nax + 2;

    double *norm1   = (double*)calloc(sz, sizeof(double));
    double *norm2   = (double*)calloc(sz, sizeof(double));
    double *norm11  = (double*)calloc(sz, sizeof(double));
    double *normn10 = (double*)calloc(sz, sizeof(double));

    /* norm1m, norm2m, normn1 are (nax+2) x (nax+2) */
    double *norm1m  = (double*)calloc(sz*sz, sizeof(double));
    double *norm2m  = (double*)calloc(sz*sz, sizeof(double));
    double *normn1  = (double*)calloc(sz*sz, sizeof(double));

    /* p is (nax+2) x (nax+4) — extra cols for safety with p(ni, ni+2) */
    int p_cols = sz + 3;
    double *p = (double*)calloc(sz * p_cols, sizeof(double));

    double *ctil = (double*)calloc(sz, sizeof(double));
    double *stil = (double*)calloc(sz, sizeof(double));

    if (!norm1 || !norm2 || !norm11 || !normn10 ||
        !norm1m || !norm2m || !normn1 || !p || !ctil || !stil)
        mexErrMsgIdAndTxt("funcGottliebNorm:memory", "Memory allocation failed.");

    /* ---- Step 1: Precompute normalisation coefficients ---- */
    /* MATLAB loop: for n = 2:nax+1 — note nax+1 not nax */
    for (int n = 2; n <= nax + 1; n++) {
        norm1[n]   = sqrt((double)(2*n+1) / (double)(2*n-1));
        norm2[n]   = sqrt((double)(2*n+1) / (double)(2*n-3));
        norm11[n]  = sqrt((double)(2*n+1) / (double)(2*n)) / (double)(2*n-1);
        normn10[n] = sqrt((double)((n+1)*n) / 2.0);
        for (int m = 1; m <= n; m++) {
            norm1m [IDX(n,m,sz)] = sqrt((double)((n-m)*(2*n+1))   / (double)((n+m)*(2*n-1)));
            norm2m [IDX(n,m,sz)] = sqrt((double)((n-m)*(n-m-1)*(2*n+1)) /
                                        (double)((n+m)*(n+m-1)*(2*n-3)));
            normn1 [IDX(n,m,sz)] = sqrt((double)((n+m+1)*(n-m)));
        }
    }

    /* ---- Step 2: Rotate position into ECEF frame x = rnp * xin ---- */
    /* rnp is 3x3, column-major */
    double x[3];
    x[0] = rnp[0]*xin_ptr[0] + rnp[3]*xin_ptr[1] + rnp[6]*xin_ptr[2];
    x[1] = rnp[1]*xin_ptr[0] + rnp[4]*xin_ptr[1] + rnp[7]*xin_ptr[2];
    x[2] = rnp[2]*xin_ptr[0] + rnp[5]*xin_ptr[1] + rnp[8]*xin_ptr[2];

    /* ---- Step 3: Derived quantities ---- */
    double r    = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double ri   = 1.0 / r;
    double xor  = x[0] * ri;
    double yor  = x[1] * ri;
    double zor  = x[2] * ri;
    double ep   = zor;
    double reor = re * ri;
    double reorn = reor;
    double muor2 = mu * ri * ri;

    /* ---- Step 4: Initialise p (Legendre polynomials) ---- */
    /* MATLAB: p(1,1)=1; p(1,2)=0; p(1,3)=0; */
    /* MATLAB 1-indexed → C 0-indexed: p(n,m) → p[IDX(n-1, m-1, sz)] */
    /* We keep 1-based indexing to match MATLAB exactly — waste row/col 0 */
    p[IDX(1,1,sz)] = 1.0;
    p[IDX(1,2,sz)] = 0.0;
    p[IDX(1,3,sz)] = 0.0;
    p[IDX(2,2,sz)] = sqrt(3.0);
    p[IDX(2,3,sz)] = 0.0;
    p[IDX(2,4,sz)] = 0.0;

    /* MATLAB: for n=2:nax; ni=n+1; p(ni,ni)=norm11(n)*p(n,n)*(2n-1); */
    for (int n = 2; n <= nax; n++) {
        int ni = n + 1;
        p[IDX(ni, ni,   sz)] = norm11[n] * p[IDX(n, n, sz)] * (double)(2*n-1);
        p[IDX(ni, ni+1, sz)] = 0.0;
        p[IDX(ni, ni+2, sz)] = 0.0;
    }

    /* ---- Step 5: Initialise ctil, stil ---- */
    ctil[1] = 1.0;
    stil[1] = 0.0;
    ctil[2] = xor;
    stil[2] = yor;

    /* ---- Step 6: Initialise sums ---- */
    double sumh  = 0.0;
    double sumgm = 1.0;
    double sumj  = 0.0;
    double sumk  = 0.0;

    /* ---- Step 7: p(2,1) ---- */
    p[IDX(2,1,sz)] = sqrt(3.0) * ep;

    /* ---- Step 8: Main loop ---- */
    for (int n = 2; n <= nax; n++) {
        int ni   = n + 1;
        int nm1  = n - 1;
        int np1  = n + 1;
        int n2m1 = 2*n - 1;

        reorn *= reor;

        /* Recurrence for p(ni, n) and p(ni, 1), p(ni, 2) */
        p[IDX(ni, n,   sz)] = normn1[IDX(n, n-1, sz)] * ep * p[IDX(ni, ni, sz)];
        p[IDX(ni, 1,   sz)] = ((double)n2m1 * ep * norm1[n]         * p[IDX(n,   1, sz)]
                               - (double)nm1       * norm2[n]         * p[IDX(nm1, 1, sz)])
                              / (double)n;
        p[IDX(ni, 2,   sz)] = ((double)n2m1 * ep * norm1m[IDX(n,1,sz)] * p[IDX(n,   2, sz)]
                               - (double)n         * norm2m[IDX(n,1,sz)] * p[IDX(nm1, 2, sz)])
                              / (double)nm1;

        /* C and S are (nax+1)x(nax+1) in MATLAB, 1-indexed */
        /* In C (column-major): C(ni, 1) = C_mat[IDX(ni-1, 0, CS_rows)] */
        double sumhn  = normn10[n] * p[IDX(ni, 2, sz)] * C_mat[IDX(ni-1, 0, CS_rows)];
        double sumgmn = p[IDX(ni, 1, sz)] * C_mat[IDX(ni-1, 0, CS_rows)] * (double)np1;

        if (mmax > 0) {
            /* Inner recurrence for p(ni, mi), m=2..n-2 */
            for (int m = 2; m <= n-2; m++) {
                int mi  = m + 1;
                p[IDX(ni, mi, sz)] = ((double)n2m1 * ep * norm1m[IDX(n,m,sz)] * p[IDX(n,   mi, sz)]
                                      - (double)(nm1+m) * norm2m[IDX(n,m,sz)] * p[IDX(nm1, mi, sz)])
                                     / (double)(n-m);
            }

            double sumjn = 0.0;
            double sumkn = 0.0;

            ctil[ni] = ctil[2]*ctil[ni-1] - stil[2]*stil[ni-1];
            stil[ni] = stil[2]*ctil[ni-1] + ctil[2]*stil[ni-1];

            int lim = (n < mmax) ? n : mmax;

            for (int m = 1; m <= lim; m++) {
                int mi  = m + 1;
                int mm1 = mi - 1;
                int mp1 = mi + 1;

                double mxpnm  = (double)m * p[IDX(ni, mi, sz)];
                double bnmtil = C_mat[IDX(ni-1, mi-1, CS_rows)] * ctil[mi]
                              + S_mat[IDX(ni-1, mi-1, CS_rows)] * stil[mi];

                sumhn  += normn1[IDX(n, m, sz)] * p[IDX(ni, mp1, sz)] * bnmtil;
                sumgmn += (double)(n+m+1) * p[IDX(ni, mi, sz)] * bnmtil;

                double bnmtm1 = C_mat[IDX(ni-1, mi-1, CS_rows)] * ctil[mm1]
                              + S_mat[IDX(ni-1, mi-1, CS_rows)] * stil[mm1];
                double anmtm1 = C_mat[IDX(ni-1, mi-1, CS_rows)] * stil[mm1]
                              - S_mat[IDX(ni-1, mi-1, CS_rows)] * ctil[mm1];

                sumjn += mxpnm * bnmtm1;
                sumkn -= mxpnm * anmtm1;
            }

            sumj += reorn * sumjn;
            sumk += reorn * sumkn;
        }

        sumh  += reorn * sumhn;
        sumgm += reorn * sumgmn;
    }

    /* ---- Step 9: Compute gravity in ECEF ---- */
    double lambda = sumgm + ep * sumh;

    double g[3];
    g[0] = -muor2 * (lambda * xor - sumj);
    g[1] = -muor2 * (lambda * yor - sumk);
    g[2] = -muor2 * (lambda * zor - sumh);

    /* ---- Step 10: Rotate back to ECI: aGravity = rnp' * g ---- */
    /* rnp is column-major, so rnp' has element [i,j] = rnp[IDX(j,i,3)] */
    double aGravity[3];
    aGravity[0] = rnp[0]*g[0] + rnp[1]*g[1] + rnp[2]*g[2];
    aGravity[1] = rnp[3]*g[0] + rnp[4]*g[1] + rnp[5]*g[2];
    aGravity[2] = rnp[6]*g[0] + rnp[7]*g[1] + rnp[8]*g[2];

    /* ---- Step 11: Create output array ---- */
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
    double *out = mxGetPr(plhs[0]);
    out[0] = aGravity[0];
    out[1] = aGravity[1];
    out[2] = aGravity[2];

    /* ---- Cleanup ---- */
    free(norm1);   free(norm2);  free(norm11); free(normn10);
    free(norm1m);  free(norm2m); free(normn1);
    free(p);       free(ctil);   free(stil);
}
