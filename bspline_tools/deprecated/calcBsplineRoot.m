function [z, c0, b] = calcBsplineRoot( N )
% Calculate the roots for the direct filter coefs (z polynomial) of an Nth
% order B-spline

directFilts = { {1,1}, {1,1}, {8, [1, 6, 1]}, {6, [1,4,1]}, ...
                {384, [1, 76, 230, 76, 1]}, {120, [1, 26, 66, 26, 1]},...
                {46080, [1, 722, 10543, 23548, 10543, 722, 1]},...
                {5040, [1, 120, 1191, 2416, 1191, 120, 1]} };

df = directFilts{N +1};
if(length(df{2}) > 1)
    c0 = df{1};
    z = roots(df{2});
    z = z((length(z)/2 + 1): end);
else
    c0 = 1;
    z = [];
end

b=df{2};

end