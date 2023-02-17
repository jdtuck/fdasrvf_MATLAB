function gam = v_to_gam(vec)
    % v_to_gam Convert shooting vectors to warping functions
    % -----------------------------------------------------_

    [n,T] = size(vec);
    time = linspace(0,1,T);

    gam = zeros(n,T);
    for i=1:n
        psi = exp_map(mu,vec(i,:));
        gam_tmp = cumtrapz(time,psi.^2);
        gam(i,:) = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
    end
