function q = f_to_srvf(f,time,smooth,spl)
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
% smooth: use smoothing splines (default: true)
% spl: use b-spline computation (old behavior)
%
% Output:
% q: matrix of SRSFs

arguments
    f double
    time double
    smooth=true
    spl=false
end

binsize = mean(diff(time));
[M, N] = size(f);

if smooth
    fy = zeros(M,N);
    for ii = 1:N
        y = fit(time(:), f(:,ii),'smoothingspline');
        fy(:,ii) = differentiate(y, time);
    end
    q = fy./sqrt(abs(fy)+eps);
elseif spl
    fy = zeros(M,N);
    for ii = 1:N
        y = Bspline(f(:,ii),3);
        ydiff = diff(y);
        fy(:,ii) = ydiff(1:length(time))/binsize;
        f(:,ii) = y(1:length(time));
    end
    q = fy./sqrt(abs(fy)+eps);
else
    if size(f,2)>1
        [~,fy] = gradient(f,binsize,binsize);
    else
        fy = gradient(f,binsize);
    end
    q = fy./sqrt(abs(fy)+eps);
end

