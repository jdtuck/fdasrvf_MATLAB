function t1 = arclength(f)
% ARCLENGTH Calculate arclength
t1 = zeros(length(f),1);
t1(1) = 0;
t1(2:end) = abs(diff(f));
t1 = cumsum(t1);
end