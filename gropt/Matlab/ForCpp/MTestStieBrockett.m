function MTestStieBrockett()
    n = 12;
    p = 4;
    B = randn(n, n);
    B = B + B';
    D = (p:-1:1)';
    Xinitial = orth(randn(n, p));
    SolverParams.method = 'RBFGS';
    SolverParams.IsCheckParams = 1;
    SolverParams.DEBUG = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieBrockett(B, D, Xinitial, HasHHR, SolverParams);
end
