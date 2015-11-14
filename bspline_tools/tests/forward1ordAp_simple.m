function cp = forward1ordAp_simple( x, zi )
% forward direction 1st order allpole filter
%
% The reflected input data is used as initial conditions

b0 = 1;
a0 = [1, -zi];

k0 = round(log(0.5*sqrt(eps)) / log(abs(zi)));
k0 = min(k0,length(x)-1);

zpows = cumprod(zi*ones(1,k0));
zpows = [1,zpows];

cp = zeros(size(x));
cp(1) = sum(zpows.*x(1:(k0+1)));
for n = 2:length(x)
    cp(n) = x(n) + zi*cp(n-1);
end

end

