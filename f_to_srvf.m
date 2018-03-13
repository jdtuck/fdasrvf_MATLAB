function q = f_to_srvf(f,time)
% F_TO_SRVF Convert function to Square-Root Velocity Function
% -------------------------------------------------------------------------
% Convert to SRSF
% 
% Usage: q = f_to_srvf(f,time)
%
% This function converts functions to srsf
%
% Input:
% f: matrix of functions
% time: vector of time samples
% 
% Output:
% q: matrix of SRSFs
binsize = mean(diff(time));
[M, N] = size(f);

fy = zeros(M,N);
for ii = 1:N
    y = Bspline(f(:,ii),3);
    ydiff = diff(y);
    fy(:,ii) = ydiff(1:length(time))/binsize;
    f(:,ii) = y(1:length(time));
end
q = fy./sqrt(abs(fy)+eps);
