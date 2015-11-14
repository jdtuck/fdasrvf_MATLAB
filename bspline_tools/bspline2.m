function y = bspline2( x )
% B-spline function of 2nd order

y = zeros(size(x));
x0range = abs(x) < 1/2;
x1range = abs(x)>=1/2 & abs(x)<3/2;
x0 = x(x0range);
x1 = x(x1range);
y(x0range) = -x0.^2 + 3/4;
y(x1range) = 0.5*(abs(x1)-1.5).^2;

end

