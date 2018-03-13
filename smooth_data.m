function f = smooth_data(f, sparam)
% SMOOTH_DATA Smooth Functions
% -------------------------------------------------------------------------
% This function smooths functions using standard box filter
%
% Usage:  = smooth_data(f, sparam)
%
% Input:
% f: matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
% sparam:number of times to run box filter
%
% Ouptut:
% fo: matrix of smoothed functions
[M,N] = size(f);
for r = 1:sparam
    for i = 1:N
        f(2:(M-1),i) = (f(1:(M-2),i)+2*f(2:(M-1),i) + f(3:M,i))/4;
    end
end
