function cn = reverse1ordAp_simple( cp, zi )
% reverse direction 1st order allpole filter
%
% The reflected input data is used as initial conditions

cn = zeros(size(cp));
cn(end) = (zi/(zi*zi-1)) * (cp(end)+zi*cp(end-1));

for n = (length(cp)-1):-1:1
    cn(n) = zi*(cn(n+1) - cp(n));
end

end

