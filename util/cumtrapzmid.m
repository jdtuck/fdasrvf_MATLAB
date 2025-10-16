function fn = cumtrapzmid(x,y,c,mid)
% CUMTRAPZMID Cummulative Trapezodial Integration from mid point
a = length(x);

% case < mid
fn = zeros(1,a);
tmpx = x(mid-1:-1:1);
tmpy = y(mid-1:-1:1);
tmp = c + cumtrapz(tmpx,tmpy);
fn(1:mid-1) = tmp(end:-1:1);

% case >= mid
fn(mid:a) = c + cumtrapz(x(mid:end),y(mid:end));
