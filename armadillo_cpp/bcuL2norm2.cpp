#include "armadillo/mex_interface/armaMex.hpp"
#include "mex.h"
#include "matrix.h"
#include "pair_align_functions_expomap.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check the number of input arguments.
  if (nrhs != 2)
    mexErrMsgTxt("usage: out = bcuL2norm2(x, y)");
  
  // Check type of input.
  if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) )
    mexErrMsgTxt("Input must be of type double.");
  
  // Check if input is real.
  if ( (mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])) )
    mexErrMsgTxt("Input must be real.");
  
  // Create matrices X and Y from the first and second argument.
  vec x = armaGetPr(prhs[0]);
  vec y = armaGetPr(prhs[1]);
  
  vec out = cuL2norm2(x, y);
  
  // Create the output argument plhs[0] to return out
  plhs[0] = armaCreateMxMatrix(out.n_rows, out.n_cols);
  
  // Return the vec out as plhs[0] in Matlab/Octave
  armaSetPr(plhs[0], out);
  
  return;
  
}