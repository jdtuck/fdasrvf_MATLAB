#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "dp_grid.h"
#include "dp_nbhd.h"


/* Signature:
 * function [G, T, dist] = DynamicProgrammingQ2( Q1, T1, Q2, T2, tv1, tv2, lam, nbhd_dim )
 * function [G, T, dist] = DynamicProgrammingQ2( Q1, T1, Q2, T2, tv1, tv2, lam, nbhd_dim, pen )
 *
 * Q1, Q2    dim x n arrays of SRVF samples (column-major).  Only the first
 *           nsamps-1 columns are read; Q is piecewise constant between
 *           changepoints.
 * T1, T2    1 x nsamps changepoint parameters.
 * tv1, tv2  1 x ntv parameter values defining the DP grid.
 * lam       scalar warping penalty weight.
 * nbhd_dim  scalar >= 1, size of the DP search neighborhood.
 * pen       optional scalar penalty type: 0 = none, 1 = roughness,
 *           2 = l2gam, 3 = l2psi, 4 = geodesic.  Defaults to 1 (roughness),
 *           which is the penalty this function used before pen existed.
 *
 * G and T are the warping function and its parameterization; dist is the
 * cost of the optimal path and is only computed into an output when asked
 * for.  */
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]){
  double *Q1 = 0;
  double *T1 = 0;
  double *Q2 = 0;
  double *T2 = 0;
  double lam;
  double pen_arg;
  int pen;
  double nbhd_arg;
  size_t nbhd_dim;
  int nsamps1;
  int nsamps2;
  double *tv1 = 0;
  double *tv2 = 0;
  int *idxv1 = 0;
  int *idxv2 = 0;
  int ntv1;
  int ntv2;
  double *G = 0;
  double *T = 0;
  int Gsize;
  int dim = 0;
  double *E = 0; /* E[ntv1*j+i] = cost of best path to (tv1[i],tv2[j]) */
  int *P = 0; /* P[ntv1*j+i] = predecessor of (tv1[i],tv2[j]) along best path */
  int Galloc_size;
  double res;
  size_t nbhd_count = 0; /* Number of indexes */
  Pair *dp_nbhd = 0;
  mxArray *Garr = 0;
  mxArray *Tarr = 0;
  int i;

  /* Inputs are validated here: this MEX file is called directly from
   * optimum_reparam.m, so there is no m-file wrapper doing it for us. */
  if ( nrhs != 8 && nrhs != 9 )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput", "Eight or nine inputs required: "
      "Q1, T1, Q2, T2, tv1, tv2, lam, nbhd_dim and optionally pen." );
  }
  if ( nlhs > 3 )
  {
    mexErrMsgIdAndTxt( "dp:InvalidOutput", "Too many output arguments." );
  }

  for ( i=0; i<nrhs; ++i )
  {
    if ( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
    {
      mexErrMsgIdAndTxt( "dp:InvalidInput",
        "Input %d must be a real, non-sparse double array.", i+1 );
    }
  }

  /* T1, T2, tv1 and tv2 are indexed as row vectors via mxGetN(). */
  for ( i=1; i<6; ++i )
  {
    if ( i==2 ) continue; /* prhs[2] is Q2 */
    if ( mxGetM(prhs[i]) != 1 )
    {
      mexErrMsgIdAndTxt( "dp:InvalidInput",
        "Input %d must be a row vector.", i+1 );
    }
  }

  if ( mxGetNumberOfElements(prhs[6]) != 1 ||
       mxGetNumberOfElements(prhs[7]) != 1 ||
       ( nrhs > 8 && mxGetNumberOfElements(prhs[8]) != 1 ) )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "lam, nbhd_dim and pen must be scalars." );
  }

  /* Guard the mwSize -> int narrowing below. */
  if ( mxGetM(prhs[0]) > INT_MAX || mxGetN(prhs[1]) > INT_MAX ||
       mxGetN(prhs[3]) > INT_MAX || mxGetN(prhs[4]) > INT_MAX ||
       mxGetN(prhs[5]) > INT_MAX )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput", "Input dimensions are too large." );
  }

  Q1 = mxGetPr( prhs[0] );
  T1 = mxGetPr( prhs[1] );
  Q2 = mxGetPr( prhs[2] );
  T2 = mxGetPr( prhs[3] );
  tv1 = mxGetPr( prhs[4] );
  tv2 = mxGetPr( prhs[5] );
  lam = mxGetScalar( prhs[6] );

  nbhd_arg = mxGetScalar( prhs[7] );
  /* The !(>=1) form also rejects NaN.  The upper bound keeps the cast to
   * size_t well defined; anything near it fails to allocate anyway. */
  if ( !(nbhd_arg >= 1.0) || nbhd_arg > 65535.0 )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "nbhd_dim must be a scalar between 1 and 65535." );
  }
  nbhd_dim = (size_t)nbhd_arg;

  /* 0 = no penalty, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic */
  pen_arg = ( nrhs > 8 ) ? mxGetScalar( prhs[8] ) : (double)DP_PEN_ROUGHNESS;
  /* Written so that NaN is rejected too. */
  if ( !(pen_arg >= 0.0) || pen_arg > 4.0 || pen_arg != floor(pen_arg) )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "pen must be one of 0, 1, 2, 3 or 4." );
  }
  pen = (int)pen_arg;

  dim = (int)mxGetM( prhs[0] );
  nsamps1 = (int)mxGetN( prhs[1] ); /* = columns(T1) */
  nsamps2 = (int)mxGetN( prhs[3] ); /* = columns(T2) */
  ntv1 = (int)mxGetN( prhs[4] );
  ntv2 = (int)mxGetN( prhs[5] );

  if ( dim < 1 || (int)mxGetM( prhs[2] ) != dim )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "Q1 and Q2 must have the same (nonzero) number of rows." );
  }
  if ( nsamps1 < 2 || nsamps2 < 2 || ntv1 < 2 || ntv2 < 2 )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "T1, T2, tv1 and tv2 must each have at least two elements." );
  }
  /* dp_edge_weight() reads Q columns 0 .. nsamps-2. */
  if ( mxGetN( prhs[0] ) < (mwSize)(nsamps1-1) ||
       mxGetN( prhs[2] ) < (mwSize)(nsamps2-1) )
  {
    mexErrMsgIdAndTxt( "dp:InvalidInput",
      "Q1 and Q2 must have at least columns(T)-1 columns." );
  }

  Galloc_size = ntv1>ntv2 ? ntv1 : ntv2;

  /* Sizes are computed in size_t: ntv1*ntv2 overflows int for large grids. */
  if ( !(idxv1=(int*)mxMalloc((size_t)ntv1*sizeof(int))) )
  {
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate idxv1" );
  }
  if ( !(idxv2=(int*)mxMalloc((size_t)ntv2*sizeof(int))) )
  {
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate idxv2" );
  }
  if ( !(E=(double*)mxMalloc((size_t)ntv1*(size_t)ntv2*sizeof(double))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate E" );
  }
  if ( !(P=(int*)mxCalloc((size_t)ntv1*(size_t)ntv2,sizeof(int))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate P" );
  }
  if ( !(Garr=mxCreateDoubleMatrix(1,Galloc_size,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
  }
  if ( !(Tarr=mxCreateDoubleMatrix(1,Galloc_size,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
  }

  G = mxGetPr( Garr );
  T = mxGetPr( Tarr );

  if ( !(dp_nbhd = dp_generate_nbhd(nbhd_dim, &nbhd_count)) )
  {
    mexErrMsgIdAndTxt( "dp:AllocFailed",
      "failed to allocate the search neighborhood for nbhd_dim=%g", nbhd_arg );
  }

  /* dp_costs() needs indexes for gridpoints precomputed */
  dp_all_indexes( T1, nsamps1, tv1, ntv1, idxv1 );
  dp_all_indexes( T2, nsamps2, tv2, ntv2, idxv2 );

  /* Compute cost of best path from (0,0) to every other grid point */
  res = dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, 
    dim, tv1, idxv1, ntv1, tv2, idxv2, ntv2, E, P, lam, pen,
	nbhd_count, dp_nbhd );

  /* Reconstruct best path from (0,0) to (1,1) */
  Gsize = dp_build_gamma( P, tv1, ntv1, tv2, ntv2, G, T );
  mxSetN( Garr, Gsize );
  mxSetN( Tarr, Gsize );

  /* plhs only holds nlhs entries, so assign nothing the caller did not ask
   * for.  plhs[0] is always writable, even when nlhs is 0. */
  plhs[0] = Garr;
  if ( nlhs > 1 )
    plhs[1] = Tarr;
  else
    mxDestroyArray( Tarr );

  if ( nlhs > 2 )
  {
    if ( !(plhs[2]=mxCreateDoubleScalar(res)) )
    {
      mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleScalar failed" );
    }
  }

  free( dp_nbhd );
  mxFree( idxv1 );
  mxFree( idxv2 );
  mxFree( E );
  mxFree( P );
}
