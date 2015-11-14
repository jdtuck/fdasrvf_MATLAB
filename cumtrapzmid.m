function fn = cumtrapzmid(x,y,c)
a = length(x);
mid = round(a/2);

% case < mid
fn = zeros(1,a);
tmpx = fliplr(x(mid-1:-1:1));
tmpy = fliplr(y(mid-1:-1:1));
tmp = c + cumtrapz(tmpx,tmpy);
fn(1:mid-1) = tmp(end:-1:1);

% case >= mid
fn(mid:a) = c + cumtrapz(x(mid:end),y(mid:end));
