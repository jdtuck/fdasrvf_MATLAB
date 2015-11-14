function y=bsplineN(x,n)
% B-spline function of order n
%


% Use specialised functions for lower order B-splines,
% and the general function for order > 4
if(n==0)
    y = bspline0(x);
elseif(n==1)
    y = bspline1(x);
elseif(n==2)
    y = bspline2(x);
elseif(n==3)
    y = bspline3(x);
elseif(n==4)
    y = bspline4(x);
else
    y = bsplineGeneralN(x,n);
end

end
    
    