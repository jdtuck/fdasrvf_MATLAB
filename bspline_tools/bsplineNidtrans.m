function y = bsplineNidtrans( x, N, boundaryfun )
% bsplineNdtrans - polynomial B-spline direct transform
%   y = bsplineNidtrans(x, N) calculates signal samples from N-order spline
%   coefficient vector x
%   y = bsplineNidtrans(x, N, bfun) uses bfun as boundary function to extend
%   x beyond its original domain. Default is @mirrorbound_1

if(nargin<3)
    boundaryfun = @mirrorbound_1;
end

[b, c] = idtrans_FIR_coefs(N);

x = filterFIR( b, x, boundaryfun );

% Apply scaling
y = x*(1/c);

end

