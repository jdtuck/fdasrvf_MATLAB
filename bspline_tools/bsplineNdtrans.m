function y = bsplineNdtrans( x, N, boundaryfun )
% bsplineNdtrans - polynomial B-spline direct transform
%   y = bsplineNdtrans(x, N) calculates N-order interpolating spline coefficients for
%   input vector x
%   y = bsplineNdtrans(x, N, bfun) uses bfun as boundary function to extend
%   x beyond its original domain. Default is @mirrorbound_1

if(nargin<3)
    boundaryfun = @mirrorbound_1;
end

[b, c] = idtrans_FIR_coefs(N);

% Find the z-transform poles of the direct transform 
z = roots(b);

% The roots appear in conjugate pairs, keep upper part
if(z)
    z = z((length(z)/2 + 1): end);
end

for zn = z'
    x = filterAPsym2ord( x, zn, boundaryfun);
end

% Apply scaling
y = x*c;

end

