function h = gam_to_h(gam, smooth)
% GAM_TO_H Convert warping function to shooting vectors
% -----------------------------------------------------

arguments
    gam double
    smooth=true
end

[T,n] = size(gam);
time = linspace(0,1,T);

h = zeros(T,n);
binsize = mean(diff(time));

if smooth
    for i = 1:n
        y = fit(time', gam(:,i),'smoothingspline','SmoothingParam',.99999);
        fy = differentiate(y, time);
        idx = fy <= 0;
        fy(idx) = eps;
        h(:,i) = log(fy) - trapz(time, log(fy));
    end
else
    for i=1:n
        psi = log(gradient(gam(:,i),binsize));
        h(:,i) = psi - trapz(time, psi);
    end
end
