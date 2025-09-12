function psi = gam_to_psi(gam, smooth)
% GAM_TO_PSI Convert warping function to hilbert sphere points
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
