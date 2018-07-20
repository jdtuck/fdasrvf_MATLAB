
%% ARMIJO
num = 1000;
LS_ratios = [1/2, 1/4, 1/16, 1/64];
n = 12; p = 8;

tab = zeros(num, length(LS_ratios), 8);
for k = 1 : length(LS_ratios)
    clearvars -except num LS_ratios k tab n p
    for j = 1 : num
        j
        randn('state', j);
        B = randn(n, n);
        B = B + B';
        D = (p:-1:1)';
        X = orth(randn(n, p));
        clear('SolverParams');
        SolverParams.method = 'RBFGS';
        SolverParams.LS_ratio1 = LS_ratios(k);
        SolverParams.LS_ratio2 = 1 - LS_ratios(k);
        SolverParams.IsCheckParams = 1;
        SolverParams.DEBUG = 1;
        [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieBrockett(B, D, X, 0, SolverParams);
        % [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, SolverParams);

        tab(j, k, 1) = iter;
        tab(j, k, 2) = nf;
        tab(j, k, 3) = ng;
        tab(j, k, 4) = nR;
        tab(j, k, 5) = nV;
        tab(j, k, 6) = nVp;
        tab(j, k, 7) = nH;
        tab(j, k, 8) = ComTime;
    end
end

dlmwrite(['RBFGSARMIJO'], tab);

%% WOLFE
LS_betas = [1-1/2, 1-1/4, 1-1/16, 1-1/64];

tab = zeros(num, length(LS_betas), 8);
for k = 1 : length(LS_betas)
    clearvars -except num LS_betas k tab n p
    for j = 1 : num
        j
        randn('state', j);
        B = randn(n, n);
        B = B + B';
        D = (p:-1:1)';
        X = orth(randn(n, p));
        clear('SolverParams');
        SolverParams.method = 'RBFGS';
        SolverParams.LS_beta = LS_betas(k);
        SolverParams.LineSearch_LS = 1;
        SolverParams.IsCheckParams = 1;
        SolverParams.DEBUG = 1;
        [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieBrockett(B, D, X, 1, SolverParams);
        % [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, SolverParams);

        tab(j, k, 1) = iter;
        tab(j, k, 2) = nf;
        tab(j, k, 3) = ng;
        tab(j, k, 4) = nR;
        tab(j, k, 5) = nV;
        tab(j, k, 6) = nVp;
        tab(j, k, 7) = nH;
        tab(j, k, 8) = ComTime;
    end
end
dlmwrite(['RBFGSWOLFE'], tab);
