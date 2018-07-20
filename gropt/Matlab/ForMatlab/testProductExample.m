function output = testProductExample()
n = 5;
p = 2;
m = 6;
q = 3;
B1 = randn(n, n); B1 = B1 + B1';
D1 = sparse(diag(p : -1 : 1));
B2 = randn(n, n); B2 = B2 + B2';
D2 = sparse(diag(p : -1 : 1));
B3 = randn(m, m); B3 = B3 + B3';
D3 = sparse(diag(q : -1 : 1));

fhandle = @(x)f(x, B1, D1, B2, D2, B3, D3);
gfhandle = @(x)gf(x, B1, D1, B2, D2, B3, D3);
Hesshandle = @(x, eta)Hess(x, eta, B1, D1, B2, D2, B3, D3);

SolverParams.method = 'RSD';
SolverParams.IsCheckParams = 1;
SolverParams.DEBUG = 1;
% SolverParams.IsCheckGradHess = 1;

ManiParams.IsCheckParams = 1;
ManiParams(1).name = 'Stiefel';
ManiParams(1).numofmani = 2;
ManiParams(1).n = n;
ManiParams(1).p = p;
ManiParams(1).ParamSet = 2;
ManiParams(2).name = 'Stiefel';
ManiParams(2).numofmani = 1;
ManiParams(2).n = m;
ManiParams(2).p = q;
HasHHR = 0;

X1 = orth(randn(n, p)); X2 = orth(randn(n, p)); X3 = orth(randn(m, q));
initialX.main = [X1(:); X2(:); X3(:)];
% initialX{1} = X1;
% initialX{2} = X2;
% initialX{3} = X3;

output = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX);

% [V1, D1]=eig(B1)
% [V2, D2]=eig(B2)
% [V3, D3]=eig(B3)
% 
% Xopt1 = reshape(output.main(1 : n * p), n, p)
% Xopt2 = reshape(output.main(n * p + 1 : 2 * n * p), n, p)
% Xopt3 = reshape(output.main(2 * n * p + 1 : 2 * n * p + m * q), m, q)
end

function [output, x] = f(x, B1, D1, B2, D2, B3, D3)
n = size(B1, 1); p = size(D1, 1);
m = size(B3, 1); q = size(D3, 1);
X1 = reshape(x.main(1 : n * p), n, p);
X2 = reshape(x.main(n * p + 1 : 2 * n * p), n, p);
X3 = reshape(x.main(2 * n * p + 1 : 2 * n * p + m * q), m, q);

x.BUD1 = B1 * X1 * D1; x.BUD2 = B2 * X2 * D2; x.BUD3 = B3 * X3 * D3;
output = X1(:)' * x.BUD1(:) + X2(:)' * x.BUD2(:) + X3(:)' * x.BUD3(:);
end

function [output, x] = gf(x, B1, D1, B2, D2, B3, D3)
output.main = [x.BUD1(:); x.BUD2(:); x.BUD3(:)];
output.main = 2 * output.main;
end

function [output, x] = Hess(x, eta, B1, D1, B2, D2, B3, D3)
n = size(B1, 1); p = size(D1, 1);
m = size(B3, 1); q = size(D3, 1);
eta1 = reshape(eta.main(1 : n * p), n, p);
eta2 = reshape(eta.main(n * p + 1 : 2 * n * p), n, p);
eta3 = reshape(eta.main(2 * n * p + 1 : 2 * n * p + m * q), m, q);
xi1 = 2 * B1 * eta1 * D1;
xi2 = 2 * B2 * eta2 * D2;
xi3 = 2 * B3 * eta3 * D3;
output.main = [xi1(:); xi2(:); xi3(:)];
end
