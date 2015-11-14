function [cs,zs,bs] = calcBsplineRootTable( )

cs = cell([1,8]);
zs = cell([1,8]);
bs = cell([1,8]);

for rn = 1:8;
    [zs{rn},cs{rn},bs{rn}] = calcBsplineRoot( rn - 1 );
end

end

