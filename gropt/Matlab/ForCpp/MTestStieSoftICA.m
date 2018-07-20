function MTestStieSoftICA()
    n = 12; p = 8; N = 64;
    Cs = randn(n, n, N);
    for i = 1 : N
        Cs(:, :, i) = Cs(:, :, i) + Cs(:, :, i)';
    end
    X = orth(randn(n, p));
    SolverParams.method = 'RBFGS';
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, HasHHR, SolverParams);
end
