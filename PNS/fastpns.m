function [resmat, PNS] = fastpns(x, n_pc, itype, alpha, R, thresh)
    % Input:
    %   x     (d+1) x n data matrix where each column is a unit vector.
    %
    %   n_pc     0  'full': use all 
    %            1  'approx': approximate based on 99% explained variance
    %
    %   itype    0  'seq.test' : (default) ordinary Principal Nested Sphere
    %                               with sequential tests.
    %            1  'small'    : Principal Nested SMALL Sphere
    %            2  'great'    : Principal Nested GREAT Sphere (radius pi/2)
    %  alpha     0.05 (default) : size of Type I error allowed for each test,
    %            could be any number between 0 and 1.
    %  R         100 (default) : number of bootsrap samples to be evaluated for
    %            the sequential test.

    if nargin == 1
        n_pc = 0;
        itype = 1; 
        alpha = 0.05;
        R = 100;
        thresh = 1e-15;
    end

    if nargin == 2
        itype = 1;
        alpha = 0.05;
        R = 100;
        thresh = 1e-15;
    end

    if nargin == 3
        alpha = 0.05;
        R = 100;
        thresh = 1e-15;
    end

    if nargin == 4
        R = 100;
        thresh = 1e-15;
    end

    if nargin == 5
        thresh = 1e-15;
    end

    [n, pdim] = size(x);
    if n_pc == 0
        n_pc = min(pdim-1, n-1);
    end
    if n_pc == 1
        K = cov(x);
        [U, S, V] = svd(K);
        s = diag(S);
        cumm_coef = cumsum(s) / sum(s);
        n_pc = find(cum_coef <= 0.99, 'first');
    end

    Xs = x';
    for i = 1:n
        Xs(i, :) = Xs(i, :) / norm(Xs(i, :));
    end

    muhat = mean(Xs, 1);
    muhat = muhat / norm(muhat);

    TT = Xs;
    for i = 1:n
        TT(i,:) = Xs(i,:) - sum(Xs(i,:) .* muhat) .* muhat;
    end

    TT = TT';
    K = cov(TT);
    [U, S, V] = svd(K);
    s = diag(S);
    pcapercent = sum(s[1:n_pc].^2 / sum(s.^2));

    fprintf('"Initial PNS subsphere dimension: %d\n', (n_pc+1));
    fprintf('Percentage of variability in PNS sequence %f\n', round(pcapercent * 100))

    ans = pcscore2sphere3(n_pc, muhat, Xs, TT, U);
    Xssubsphere = ans';

    resmat, PNS = PNSmainHDLSS(Xssubsphere, itype, a, R, thresh);

    PNS.spheredata = Xssubsphere;
    PNS.pca = U;
    PNS.muhat = muhat;
    PNS.n_pc = n_pc;
