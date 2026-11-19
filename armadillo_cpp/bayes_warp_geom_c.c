/*
 * bayes_warp_geom_c.c  Compiled kernel for bayes_warp_geom.m
 *
 * [nll, grad, SSE, gam, q2w] = bayes_warp_geom_c(c, Phi, rep, q1, q2, var1)
 *
 * Negative log-likelihood and exact gradient (w.r.t. basis coefficients c)
 * of q1 = (q2 o gam) sqrt(gam') + N(0, var1) on the uniform grid
 * linspace(0,1,T). rep = 0: shooting vector (v-transform), rep = 1: clr
 * (h-transform). See bayes_warp_geom.m for the model and derivation; this
 * file implements exactly the same arithmetic.
 */
#include "mex.h"
#include <math.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgTxt("usage: [nll, grad, SSE, gam, q2w] = bayes_warp_geom_c(c, Phi, rep, q1, q2, var1)");
    for (int i = 0; i < 6; i++)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
            mexErrMsgTxt("inputs must be real double");

    const mwSize T = mxGetM(prhs[1]);
    const mwSize D = mxGetN(prhs[1]);
    if (mxGetNumberOfElements(prhs[0]) != D)
        mexErrMsgTxt("length(c) must equal size(Phi,2)");
    if (mxGetNumberOfElements(prhs[3]) != T || mxGetNumberOfElements(prhs[4]) != T)
        mexErrMsgTxt("q1 and q2 must have size(Phi,1) elements");
    if (T < 2)
        mexErrMsgTxt("need at least two grid points");

    const double *c = mxGetPr(prhs[0]);
    const double *Phi = mxGetPr(prhs[1]);
    const int rep = (int) mxGetScalar(prhs[2]);
    const double *q1 = mxGetPr(prhs[3]);
    const double *q2 = mxGetPr(prhs[4]);
    const double var1 = mxGetScalar(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *grad_mx = mxCreateDoubleMatrix(D, 1, mxREAL);
    mxArray *sse_mx = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray *gam_mx = mxCreateDoubleMatrix(T, 1, mxREAL);
    mxArray *q2w_mx = mxCreateDoubleMatrix(T, 1, mxREAL);
    double *grad = mxGetPr(grad_mx);
    double *gam = mxGetPr(gam_mx);
    double *q2w = mxGetPr(q2w_mx);

    const double dt = 1.0 / (double) (T - 1);
    double *work = (double *) mxMalloc(8 * T * sizeof(double));
    double *v = work, *psi = work + T, *apsi = work + 2*T, *sgn = work + 3*T;
    double *q2g = work + 4*T, *q2dg = work + 5*T, *res = work + 6*T, *beta = work + 7*T;

    /* v = Phi * c */
    for (mwSize k = 0; k < T; k++) v[k] = 0.0;
    for (mwSize j = 0; j < D; j++) {
        const double cj = c[j];
        const double *col = Phi + j*T;
        for (mwSize k = 0; k < T; k++) v[k] += col[k] * cj;
    }

    double r = 0.0, s = 1.0, ds = 0.0, sr = 0.0;
    if (rep == 0) {
        for (mwSize k = 0; k < T; k++) {
            double wk = (k == 0 || k == T-1) ? dt/2 : dt;
            r += wk * v[k] * v[k];
        }
        r = sqrt(r);
        if (r >= M_PI) {
            *mxGetPr(plhs[0]) = mxGetInf();
            *mxGetPr(sse_mx) = mxGetInf();
            for (mwSize k = 0; k < T; k++) { gam[k] = k*dt; q2w[k] = q2[k]; }
            mxFree(work);
            if (nlhs > 1) plhs[1] = grad_mx; else mxDestroyArray(grad_mx);
            if (nlhs > 2) plhs[2] = sse_mx; else mxDestroyArray(sse_mx);
            if (nlhs > 3) plhs[3] = gam_mx; else mxDestroyArray(gam_mx);
            if (nlhs > 4) plhs[4] = q2w_mx; else mxDestroyArray(q2w_mx);
            return;
        }
        double cr = 1.0;
        if (r >= 1e-10) {
            s = sin(r) / r;
            ds = (r*cos(r) - sin(r)) / (r*r);
            sr = sin(r);
            cr = cos(r);
        }
        for (mwSize k = 0; k < T; k++) {
            psi[k] = cr + s * v[k];
            sgn[k] = (psi[k] < 0) ? -1.0 : 1.0;
            apsi[k] = fabs(psi[k]);
        }
    } else {
        double vmax = v[0];
        for (mwSize k = 1; k < T; k++) if (v[k] > vmax) vmax = v[k];
        for (mwSize k = 0; k < T; k++) {
            psi[k] = exp((v[k] - vmax) / 2);
            sgn[k] = 1.0;
            apsi[k] = psi[k];
        }
    }

    /* gam = cumtrapz(psi^2) / total */
    gam[0] = 0.0;
    for (mwSize k = 1; k < T; k++)
        gam[k] = gam[k-1] + dt/2 * (psi[k-1]*psi[k-1] + psi[k]*psi[k]);
    const double GT = gam[T-1];
    for (mwSize k = 0; k < T; k++) gam[k] /= GT;
    gam[T-1] = 1.0;

    /* interpolate q2 at gam; slope of the linear interpolant */
    const double a = 1.0 / sqrt(GT);
    double SSE = 0.0, rq = 0.0;
    for (mwSize k = 0; k < T; k++) {
        double pos = gam[k] * (T - 1);
        if (pos < 0) pos = 0;
        if (pos > T - 1) pos = T - 1;
        mwSize idx = (mwSize) floor(pos);
        if (idx > T - 2) idx = T - 2;
        double fr = pos - idx;
        q2g[k] = q2[idx]*(1 - fr) + q2[idx+1]*fr;
        q2dg[k] = (q2[idx+1] - q2[idx]) * (T - 1);
        q2w[k] = q2g[k] * apsi[k] * a;
        res[k] = q1[k] - q2w[k];
        SSE += res[k]*res[k];
        rq += res[k]*q2w[k];
    }
    *mxGetPr(plhs[0]) = SSE / (2*var1) + 0.5 * (double) T * log(var1);
    *mxGetPr(sse_mx) = SSE;

    if (nlhs > 1) {
        const double kappa = rq / (2*GT);
        double Agam = 0.0;
        for (mwSize k = 0; k < T; k++) {
            double A = res[k]*a * q2dg[k] * apsi[k];
            Agam += A * gam[k];
            q2dg[k] = A;               /* reuse: q2dg now holds A */
        }
        /* suffix sums of A, then weights on d(psi^2) */
        double Ssuf = 0.0;               /* S_{k+1} while scanning backward */
        for (mwSize kk = T; kk-- > 0; ) {
            double Sk1 = Ssuf;             /* S_{k+1} (0 for k = T-1) */
            Ssuf += q2dg[kk];              /* S_k */
            double coefP = 0.0;
            if (kk < T - 1) coefP += dt/2 * Sk1;
            if (kk > 0) coefP += dt/2 * Ssuf;
            double wk = (kk == 0 || kk == T-1) ? dt/2 : dt;
            double alpha = coefP / GT - (Agam / GT + kappa) * wk;
            double Bv = res[kk]*a * q2g[kk] * sgn[kk];
            beta[kk] = 2*alpha*psi[kk] + Bv;
        }
        if (rep == 0) {
            double sb = 0.0, sbv = 0.0;
            for (mwSize k = 0; k < T; k++) { sb += beta[k]; sbv += beta[k]*v[k]; }
            double fac = (r >= 1e-10) ? (-sr*sb + ds*sbv) / r : 0.0;
            for (mwSize j = 0; j < D; j++) {
                const double *col = Phi + j*T;
                double pb = 0.0, pwv = 0.0;
                for (mwSize k = 0; k < T; k++) {
                    double wk = (k == 0 || k == T-1) ? dt/2 : dt;
                    pb += col[k]*beta[k];
                    pwv += col[k]*wk*v[k];
                }
                grad[j] = -(s*pb + fac*pwv) / var1;
            }
        } else {
            for (mwSize j = 0; j < D; j++) {
                const double *col = Phi + j*T;
                double pb = 0.0;
                for (mwSize k = 0; k < T; k++) pb += col[k]*beta[k]*psi[k];
                grad[j] = -0.5*pb / var1;
            }
        }
    }

    mxFree(work);
    if (nlhs > 1) plhs[1] = grad_mx; else mxDestroyArray(grad_mx);
    if (nlhs > 2) plhs[2] = sse_mx; else mxDestroyArray(sse_mx);
    if (nlhs > 3) plhs[3] = gam_mx; else mxDestroyArray(gam_mx);
    if (nlhs > 4) plhs[4] = q2w_mx; else mxDestroyArray(q2w_mx);
}
