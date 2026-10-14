function gam = h_to_gam(h, parallel, nworkers)
% h_to_gam Convert shooting vectors to warping functions
% -----------------------------------------------------_

arguments
    h double
    parallel=false
    nworkers=0
end

if size(h,2) > 1
    [T, n] = size(h);
    time = linspace(0,1,T);

    gam = zeros(T, n);
    M = parfor_workers(parallel, nworkers);
    parfor (i=1:n, M)
        gam_tmp = cumtrapz(time,exp(h(:,i)));
        gam(:,i) = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
    end
else
    T = length(h);
    time = linspace(0,1,T);
    gam_tmp = cumtrapz(time,exp(h));
    gam = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
end
