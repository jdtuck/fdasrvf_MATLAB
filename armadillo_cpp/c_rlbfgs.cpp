#include "armadillo/mex_interface/armaMex.hpp"
#include "mex.h"
#include "matrix.h"
#include "rbfgs.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check the number of input arguments.
  if (nrhs != 6)
    mexErrMsgTxt("usage: gam = c_rlbfgs(q1, q2, time, maxiter, lam, penalty)");
  
  // Check type of input.
  if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) )
    mexErrMsgTxt("Input must be of type double.");
  
  // Check if input is real.
  if ( (mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])) )
    mexErrMsgTxt("Input must be real.");
  
  // Create matrices arguments.
  vec q1 = armaGetPr(prhs[0]);
  vec q2 = armaGetPr(prhs[1]);
  vec time = armaGetPr(prhs[2]);
  int maxiter = mxGetScalar(prhs[3]);
  double lam = mxGetScalar(prhs[4]);
  int penalty = mxGetScalar(prhs[5]);
  
  vec gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty);
  
  // Create the output argument plhs[0] to return out
  plhs[0] = armaCreateMxMatrix(gam.n_rows, gam.n_cols);
  
  // Return the vec gam as plhs[0] in Matlab/Octave
  armaSetPr(plhs[0], gam);
  
  return;
  
}