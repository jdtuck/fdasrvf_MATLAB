function y = bspline0( x )
% B-spline function of 0th order

y = zeros(size(x));
y(x>=-0.5 & x<0.5) = 1;


end

