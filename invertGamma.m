function gamI = invertGamma(gam)

N = length(gam);
x = (0:N-1)/(N-1);
gamI = interp1(gam,x,x);
if isnan(gamI(N))
    gamI(N) = 1;
else
    gamI = gamI./gamI(N);
end

