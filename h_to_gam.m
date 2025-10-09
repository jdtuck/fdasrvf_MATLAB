function gam = h_to_gam(h)
    % h_to_gam Convert shooting vectors to warping functions
    % -----------------------------------------------------_

    if size(h,2) > 1
        [n,T] = size(h);
        time = linspace(0,1,T);
    
        gam = zeros(n,T);
        mu = ones(1, T);
        for i=1:n
            gam_tmp = cumtrapz(time,exp(h(i,:)));
            gam(i,:) = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
        end
    else
        T = length(h);
        time = linspace(0,1,T);
        
        mu = ones(1,T);
        gam_tmp = cumtrapz(time,exp(h));
        gam = (gam_tmp-min(gam_tmp))/(max(gam_tmp)-min(gam_tmp));
    end
