function q = f_to_srvf(f,time)
% Convert to SRSF
%
% This function converts functions to srsf
%
% @param f matrix of functions
% @param time time
% @return q matrix of SRSFs
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
