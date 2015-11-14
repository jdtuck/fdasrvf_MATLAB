function c = bspline3dtrans( x )
% Direct transform specialised for 3rd order bspline

z1 = -2 + sqrt(3);
c1 = 6;

c = c1*filterAPsym2ord( x, z1);

end

