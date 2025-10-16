function tol = mvtol_region(x, alpha, P, B )
% MVTOL_REGION Computes tolerance factor for multivariate normal
% -------------------------------------------------------------------------
% Krishnamoorthy, K. and Mondal, S. (2006), Improved Tolerance Factors for Multivariate Normal
% Distributions, Communications in Statistics - Simulation and Computation, 35, 461–478.
%
% Usage: tol = mvtol_region(x, alpha, P, B )
%
% Inputs:
% x (M,N): matrix defining N variables of M samples
% alpha: confidence level
% P: coverage level
% B: number of bootstrap samples
%
% Outputs:
% tol: tolerance factor
[n,p] = size(x);

q_squared = chi2rnd(1, [p, B])./n;
L = arrayfun(@(x){eig(rwishart(n-1,p))},1:B);
L = cell2mat(L);
c1 = sum((1+q_squared)./L);
c2 = sum((1+2.*q_squared)./(L.^2));
c3 = sum((1+3.*q_squared)./(L.^3));
a = (c2.^3)./(c3.^2);
T = (n-1)*(sqrt(c2./a) .* (chi2inv(P, a) - a) + c1);
tol = quantile(T,1-alpha);
end

function R = rwishart(df, p)
R = zeros(p,p);
R(1:p+1:end) = sqrt(chi2rnd(df:-1:(df-p+1),1,p));
if (p > 1)
    pseq = 1:(p-1);
    R(repelem(p*pseq,pseq)+cell2mat(arrayfun(@(x){1:x},pseq))) = randn(p*(p-1)/2,1);
end
R = R.' * R;
end
