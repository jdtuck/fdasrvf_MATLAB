function gamI = invertGamma(gam)
% Invert Warping Function
%
% This function calculates the inverse of gamma
%
% @param gam vector of \eqn{N} samples
% @return Returns gamI inverted vector

N = length(gam);
x = (0:N-1)/(N-1);
gamI = interp1(gam,x,x);
if isnan(gamI(N))
    gamI(N) = 1;
else
    gamI = gamI./gamI(N);
end
