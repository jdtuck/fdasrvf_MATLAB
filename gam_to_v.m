function vec = gam_to_v(gam, smooth)
% GAM_TO_V Convert warping function to shooting vectors
% -----------------------------------------------------

arguments
    gam double
    smooth=true
end

[T,n] = size(gam);
time = linspace(0,1,T);

psi = zeros(T,n);
binsize = mean(diff(time));

if smooth
    for i = 1:n
        y = fit(time', gam(:,i),'smoothingspline','SmoothingParam',.9999);
        fy = differentiate(y, time);
        idx = fy <= 0;
        fy(idx) = 0;
        psi(:,i) = sqrt(fy);
    end
else
    for i=1:n
        psi(:,i) = sqrt(gradient(gam(:,i),binsize));
    end
end

mu = ones(1, T);
vec = zeros(T,n);
for i = 1:n
    vec(:,i) = inv_exp_map(mu,psi(:,i));
end
