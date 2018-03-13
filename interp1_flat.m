function yy = interp1_flat(x,y,xx)
% INTERP1_FLAT Flat linear interpolation
flat = find(diff(x)==0);
n = length(flat);

if n==0
    yy = interp1(x,y,xx);
else
    
    yy = zeros(size(xx));
    
    i1 = 1;
    if flat(1)==1
        i2 = 1;
        j = xx==x(i2);
        yy(j) = min(y(i2:i2+1));
    else
        i2 = flat(1);
        j = xx>=x(i1) & xx<=x(i2);
        yy(j) = interp1(x(i1:i2),y(i1:i2),xx(j));
        i1 = i2;
    end
    for k=2:n
        i2 = flat(k);
        if i2>i1+1
            j = xx>=x(i1) & xx<=x(i2);
            yy(j) = interp1(x(i1+1:i2),y(i1+1:i2),xx(j));
        end
        j = xx==x(i2);
        yy(j) = min(y(i2:i2+1));
        i1 = i2;
    end
    i2 = length(x);
    j = xx>=x(i1) & xx<=x(i2);
    if i1+1==i2
        yy(j) = y(i2);
    else
        yy(j) = interp1(x(i1+1:i2),y(i1+1:i2),xx(j));
    end
end
end