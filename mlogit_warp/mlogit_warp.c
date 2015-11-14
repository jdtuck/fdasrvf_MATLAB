#include "mex.h"
#include <math.h>
#include "mlogit_warp_grad.h"

/* Signature:
 * function gamma = mlogit_warp(alpha, beta, time, q, y, gami, max_iter, tol, delta, display)
 * TO DO: Argument Checking  */

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]){
    int m1, m2;
    int *y, *max_itri, *displayi;
    double *alpha, * beta, *ti, *gami, *q, *toli, *deltai, *gamout;
    int Galloc_size;

    alpha = mxGetPr( prhs[0] );
    beta = mxGetPr( prhs[1] );
    ti = mxGetPr( prhs[2] );
    q = mxGetPr( prhs[3] );
    y = mxGetData( prhs[4] );
    gami = mxGetPr( prhs[5] );
    max_itri = mxGetData( prhs[6] );
    toli = mxGetPr( prhs[7] );
    deltai = mxGetPr( prhs[8] );
    displayi = mxGetData( prhs[9] );

    m1 = mxGetM(prhs[2]);
    if (m1 == 1){
        m1 = mxGetN(prhs[2]);
    }
    m2 = mxGetN(prhs[4]);
    Galloc_size = m1;

    if ( !(plhs[0]=mxCreateDoubleMatrix(Galloc_size,1,mxREAL)) )
    {
        mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
    }

    gamout = mxGetPr( plhs[0] );

    mlogit_warp_grad(&m1, &m2, alpha, beta, ti, gami, q, y, max_itri, toli, deltai, displayi, gamout);

}
