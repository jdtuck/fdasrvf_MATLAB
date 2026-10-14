function h = gam_to_h(gam, smooth, parallel, nworkers)
% GAM_TO_H Convert warping function to shooting vectors
% -----------------------------------------------------

arguments
    gam double
    smooth=true
    parallel=false
    nworkers=0
end

[T,n] = size(gam);
time = linspace(0,1,T);

h = zeros(T,n);
binsize = mean(diff(time));

M = parfor_workers(parallel, nworkers);
if smooth
    parfor (i = 1:n, M)
        y = fit(time', gam(:,i),'smoothingspline','SmoothingParam',.99999);
        fy = differentiate(y, time);
        idx = fy <= 0;
        fy(idx) = eps;
        h(:,i) = log(fy) - trapz(time, log(fy));
    end
else
    parfor (i=1:n, M)
        psi = log(gradient(gam(:,i),binsize));
        h(:,i) = psi - trapz(time, psi);
    end
end
