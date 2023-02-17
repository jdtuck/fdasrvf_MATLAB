function vec = gam_to_v(gam)
    % GAM_TO_V Convert warping function to shooting vectors
    % -----------------------------------------------------

    [n,T] = size(gam);
    time = linspace(0,1,T);

    psi = zeros(n,T);
    binsize = mean(diff(time));
    for i=1:n
        psi(i,:) = sqrt(gradient(gam(i,:),binsize));
    end

    mu = ones(1, TT);
    vec = zeros(n,T);
    for i = 1:n
        vec(i,:) = inv_exp_map(mu,psi(i,:));
    end
