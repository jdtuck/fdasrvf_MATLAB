function [ext, d] = extrema_1s(t,q)
% EXTREMA_1s Finds extrema in SIMUL Warping
q = round(q);

if q(1)~=0
    d = -q(1);
else
    d = q(q~=0);
    d = d(1);
end

ext = find(diff(q))+1;

ext(2:end+1) = ext;
ext(1) = 1;
ext(end+1) = length(t);
end
