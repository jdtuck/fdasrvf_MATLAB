function y = bspline3( x )
% B-spline function of 3rd order

y = zeros(size(x));
x0range = abs(x) < 1;
x1range = abs(x)>=1 & abs(x)<2;
x0 = x(x0range);
x1 = x(x1range);
y(x0range) = 2/3 - abs(x0).^2 + (1/2)*abs(x0).^3;
y(x1range) = (1/6)*(2-abs(x1)).*(2-abs(x1)).*(2-abs(x1));

end

