function approx1 = fastPNSe2s(res, PNS)
    % Fast PNS coordinate transform from Euclidean-type residual matrix to Sphere

    % Usage: Spheredata = fastPNSe2s(res, PNS)

    %  where 'res' is d x m data matrix in PNS coordinate system (for any
    %  m >= 1), 'PNS' is the structural array

    GG = PNSe2s(res, PNS);
    n = size(GG, 2);
    muhat = PNS.muhat;
    n_pc = PNS.n_pc;
    s = acos(GG(1,:));
    ones_vec = ones(n,1);
    approx1 = GG(2:(n_pc+1), :)' * PNS.pca(:, 1:n_pc)' + diag(cos(s)) * ones * muhat' / norm(muhat);
