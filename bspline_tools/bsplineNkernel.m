function b = bsplineNkernel(k,n,m)
% centered B-spline kernel of order n, expansion factor m

b = bsplineN(k/m,n);

end