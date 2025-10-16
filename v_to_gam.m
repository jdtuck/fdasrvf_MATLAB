function gam = v_to_gam(vec)
% v_to_gam Convert shooting vectors to warping functions
% -----------------------------------------------------_

if size(vec,2) > 1
    [T, n] = size(vec);
    time = linspace(0,1,T);

    gam = zeros(T,n);
    mu = ones(1, T);
    for i=1:n
        psi = exp_map(mu,vec(:,i));
        gam_tmp = cumtrapz(time,psi.^2);
        gam(:,i) = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
    end
else
    T = length(vec);
    time = linspace(0,1,T);

    mu = ones(1,T);
    psi = exp_map(mu,vec);
    gam_tmp = cumtrapz(time,psi.^2);
    gam = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
end
