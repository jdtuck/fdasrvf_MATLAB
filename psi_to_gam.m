function gam = psi_to_gam(psi)
    % psi_to_gam Convert shooting vectors to warping functions
    % -----------------------------------------------------_

    if size(psi,2) > 1
        [T, n] = size(psi);
        time = linspace(0,1,T);
    
        gam = zeros(T,n);
        for i=1:n
            gam_tmp = cumtrapz(time,psi(:,i).^2);
            gam(:,i) = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
        end
    else
        T = length(psi);
        time = linspace(0,1,T);
        
        gam_tmp = cumtrapz(time,psi.^2);
        gam = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
    end
