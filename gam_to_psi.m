function psi = gam_to_psi(gam, smooth, parallel, nworkers)
% GAM_TO_PSI Convert warping function to hilbert sphere points
% -----------------------------------------------------
% psi = sqrt(gamdot), rescaled to unit L2 norm so it lies on the Hilbert
% sphere (smoothing and clipping negative slopes otherwise move it off)

arguments
    gam double
    smooth=true
    parallel=false
    nworkers=0
end

[T,n] = size(gam);
time = linspace(0,1,T);

psi = zeros(T,n);
binsize = mean(diff(time));

M = parfor_workers(parallel, nworkers);
if smooth
    parfor (i = 1:n, M)
        y = fit(time', gam(:,i),'smoothingspline','SmoothingParam',.9999);
        fy = differentiate(y, time);
        idx = fy <= 0;
        fy(idx) = 0;
        psi(:,i) = sqrt(fy);
    end
else
    parfor (i=1:n, M)
        psi(:,i) = sqrt(max(gradient(gam(:,i),binsize),0));
    end
end

nrm = sqrt(trapz(time, psi.^2, 1));
nrm(nrm == 0) = 1;
psi = psi ./ nrm;
