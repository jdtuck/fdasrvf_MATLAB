function gamI = invertGamma(gam)
% INVERTGAMMA Invert Warping Function
% -------------------------------------------------------------------------
% This function calculates the inverse of gamma
% 
% Usage: gamI = invertGamma(gam)
%
% Input:
% gam: vector of \eqn{N} samples
%
% Output:
% gamI: inverted vector

N = length(gam);
x = (0:N-1)/(N-1);
gamI = interp1(gam,x,x);
if isnan(gamI(N))
    gamI(N) = 1;
else
    gamI = gamI./gamI(N);
end
