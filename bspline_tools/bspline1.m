function y = bspline1( x )
% B-spline function of 1st order

y = zeros(size(x));

x0range = abs(x) < 1;
x0 = x(x0range);
y(x0range) = 1-abs(x0);


end

