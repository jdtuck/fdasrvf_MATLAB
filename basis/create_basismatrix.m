function basismat = create_basismatrix(x, nbasis, norder)
% CREATE_BASISMATRIX Create Basis Matrix for BSPLINE basis
% -------------------------------------------------------------------------
% Evaluates a B-spline basis of NBASIS functions of order NORDER on
% LENGTH(X) equally spaced points spanning the range of X. The knot sequence
% is NBASIS-NORDER+2 equally spaced breakpoints with NORDER-fold
% multiplicity at each end of the range.
%
% Usage: basismat = create_basismatrix(x, nbasis, norder)
%
% Input:
% x: vector of argument values; only its range and length are used
% nbasis: number of basis functions
% norder: order of the B-splines, i.e. degree + 1 (default: 4, cubic)
%
% Output:
% basismat: (length(x) x nbasis) matrix of basis function values

arguments
    x double
    nbasis (1,1) double
    norder (1,1) double = 4
end

nbreaks = nbasis - norder + 2;
if nbreaks < 2
    error('fdasrvf:create_basismatrix:tooFewBasis', ...
        'nbasis (%d) must be at least norder (%d).', nbasis, norder)
end

rangex = [min(x) max(x)];
knots = augknt(linspace(rangex(1), rangex(2), nbreaks), norder);
basismat = spcol(knots, norder, linspace(rangex(1), rangex(2), numel(x)).');
end
