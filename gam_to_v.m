function vec = gam_to_v(gam, smooth, parallel, nworkers)
% GAM_TO_V Convert warping function to shooting vectors
% -----------------------------------------------------

arguments
    gam double
    smooth=false
    parallel=false
    nworkers=0
end

[T,n] = size(gam);

psi = gam_to_psi(gam, smooth, parallel, nworkers);

mu = ones(1, T);
vec = zeros(T,n);
M = parfor_workers(parallel, nworkers);
parfor (i = 1:n, M)
    vec(:,i) = inv_exp_map(mu,psi(:,i));
end
