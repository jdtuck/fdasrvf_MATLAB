function y = bspline4( x )
% B-spline function of 4th order

y = zeros(size(x));
x0range = abs(x) < 0.5;
x1range = abs(x)>=0.5 & abs(x)<1.5;
x2range = abs(x)>=1.5 & abs(x)<2.5;
x0 = x(x0range);
x1 = x(x1range);
x2 = x(x2range);
y(x0range) = (1/24)*(6*x0.^4 - 15*x0.^2 + 14.375);
y(x1range) = (1/24)*(-4*x1.^4 + 20*abs(x1).^3 - 30*x1.^2 + 5*abs(x1) + 13.75);
y(x2range) = (1/24)*(x2.^4 - 10*abs(x2).^3 + 37.5*x2.^2 - 62.5*abs(x2) + 39.0625);
end

