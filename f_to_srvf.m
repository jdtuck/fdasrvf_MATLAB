function q = f_to_srvf(f,time,smooth,spl,parallel,nworkers)
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
% smooth: use smoothing splines instead of the interpolating spline (default: false)
% spl: deprecated and ignored (kept for backward compatibility); when
%      smooth is false the derivative always comes from an interpolating
%      cubic spline
% paralell: compute in parallel
% nworkers: number of workers for the parallel pool (0: prompt)
%
% Output:
% q: matrix of SRSFs
%
% Note: with smooth=false the derivative is that of the interpolating cubic
% spline, which srvf_to_f integrates exactly, so f -> q -> f is accurate to
% O(h^4) for smooth f. smooth=true applies a smoothing spline, so
% high-frequency content is removed on purpose and is not recovered.

arguments
    f double
    time double
    smooth=false
    spl=false
    parallel=false
    nworkers=0
end

if parallel == 1
    start_parpool(nworkers);
end

[M, N] = size(f);
fy = zeros(M,N);

if smooth
    if parallel
        parfor ii = 1:N
            y = fit(time(:), f(:,ii),'smoothingspline');
            fy(:,ii) = differentiate(y, time);
        end
    else
        for ii = 1:N
            y = fit(time(:), f(:,ii),'smoothingspline');
            fy(:,ii) = differentiate(y, time);
        end
    end
else
    tt = time(:);
    if parallel
        parfor ii = 1:N
            fy(:,ii) = fnval(fnder(csapi(tt, f(:,ii))), tt);
        end
    else
        for ii = 1:N
            fy(:,ii) = fnval(fnder(csapi(tt, f(:,ii))), tt);
        end
    end
end
q = fy./sqrt(abs(fy)+eps);
