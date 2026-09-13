#include "armadillo/mex_interface/armaMex.hpp"
#include <mex.h>
#include <matrix.h>
#include <climits>
#include <cmath>
#include <exception>
#include "rbfgs.hpp"

namespace {

// Check an argument that must hold a real, dense, double vector.
void checkVectorArg(const mxArray *arr, const char *name) {
  if (mxGetClassID(arr) != mxDOUBLE_CLASS)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidType",
                      "%s must be of type double.", name);

  if (mxIsComplex(arr))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidType",
                      "%s must be real.", name);

  if (mxIsSparse(arr))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidType",
                      "%s must be dense.", name);

  if (mxGetNumberOfDimensions(arr) != 2 ||
      (mxGetM(arr) != 1 && mxGetN(arr) != 1))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidSize",
                      "%s must be a vector.", name);

  if (mxGetNumberOfElements(arr) < 2)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidSize",
                      "%s must have at least two elements.", name);
}

// Read an argument that must hold a real, finite, numeric scalar.
double getScalarArg(const mxArray *arr, const char *name) {
  if ((!mxIsNumeric(arr) && !mxIsLogical(arr)) || mxIsComplex(arr) ||
      mxIsSparse(arr) || mxGetNumberOfElements(arr) != 1)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidType",
                      "%s must be a real numeric scalar.", name);

  const double val = mxGetScalar(arr);

  if (!mxIsFinite(val))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidValue",
                      "%s must be finite.", name);

  return val;
}

// Copy a MATLAB vector into an Armadillo column vector, row or column alike.
vec getVectorArg(const mxArray *arr) {
  const mat tmp = armaGetPr(arr);

  return vec(tmp.memptr(), tmp.n_elem);
}

}  // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check the number of input and output arguments.
  if (nrhs != 6)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidNumInputs",
                      "usage: gam = c_rlbfgs(q1, q2, time, maxiter, lam, penalty)");

  if (nlhs > 1)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidNumOutputs",
                      "Too many output arguments.");

  // Check type and shape of the vector inputs.
  checkVectorArg(prhs[0], "q1");
  checkVectorArg(prhs[1], "q2");
  checkVectorArg(prhs[2], "time");

  const mwSize n = mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1]) != n || mxGetNumberOfElements(prhs[2]) != n)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:sizeMismatch",
                      "q1, q2 and time must have the same number of elements.");

  // Check type and range of the scalar inputs.
  const double maxiter_in = getScalarArg(prhs[3], "maxiter");
  const double lam = getScalarArg(prhs[4], "lam");
  const double penalty_in = getScalarArg(prhs[5], "penalty");

  if (maxiter_in < 1 || maxiter_in > INT_MAX ||
      maxiter_in != std::floor(maxiter_in))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidValue",
                      "maxiter must be a positive integer.");

  if (lam < 0)
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidValue",
                      "lam must be non-negative.");

  if (penalty_in < 0 || penalty_in > 3 || penalty_in != std::floor(penalty_in))
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:invalidValue",
                      "penalty must be 0 (roughness), 1 (l2gam), 2 (l2psi) or "
                      "3 (geodesic).");

  const int maxiter = static_cast<int>(maxiter_in);
  const int penalty = static_cast<int>(penalty_in);

  // Armadillo reports conformance failures by throwing, and an exception must
  // not escape mexFunction, so report them as MATLAB errors instead.
  vec gam;
  try {
    const vec q1 = getVectorArg(prhs[0]);
    const vec q2 = getVectorArg(prhs[1]);
    const vec time = getVectorArg(prhs[2]);

    gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:solverFailed",
                      "rlbfgs solver failed: %s", e.what());
  } catch (...) {
    mexErrMsgIdAndTxt("fdasrvf:c_rlbfgs:solverFailed",
                      "rlbfgs solver failed: unknown exception.");
  }

  // Create the output argument plhs[0] to return out
  plhs[0] = armaCreateMxMatrix(gam.n_rows, gam.n_cols);

  // Return the vec gam as plhs[0] in Matlab/Octave
  armaSetPr(plhs[0], gam);

  return;

}
