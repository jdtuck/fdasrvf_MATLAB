function [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = testSimpleExample()
n = 5;
p = 2;
B = randn(n, n); B = B + B';
D = sparse(diag(p : -1 : 1));
fhandle = @(x)f(x, B, D);
gfhandle = @(x)gf(x, B, D);
Hesshandle = @(x, eta)Hess(x, eta, B, D);

SolverParams.method = 'RSD';
SolverParams.IsCheckParams = 1;
SolverParams.DEBUG = 1;
% SolverParams.IsCheckGradHess = 1;

ManiParams.IsCheckParams = 1;
ManiParams.name = 'Stiefel';
ManiParams.n = n;
ManiParams.p = p;
ManiParams.ParamSet = 1;
HasHHR = 0;

initialX.main = orth(randn(n, p));

[FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX);
end

function [output, x] = f(x, B, D)
x.BUD = B * x.main * D;
output = x.main(:)' * x.BUD(:);
end

function [output, x] = gf(x, B, D)
output.main = 2 * x.BUD;
end

function [output, x] = Hess(x, eta, B, D)
output.main = 2 * B * eta.main * D;
end
