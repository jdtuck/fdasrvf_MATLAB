function [gam, psi] = bayes_coef_to_gam(c, Phi, rep)
% BAYES_COEF_TO_GAM Warping functions from basis coefficients
% -------------------------------------------------------------------------
% Usage: [gam, psi] = bayes_coef_to_gam(c, Phi, rep)
%
% Input:
% c: (D,K) coefficients, one column per warping function
% Phi: (T,D) basis on the uniform grid linspace(0,1,T)
% rep: 0 = shooting vector (v-transform), 1 = clr (h-transform)
%
% Output:
% gam: (T,K) warping functions
% psi: (T,K) SRSFs sqrt(gam') of the warping functions (unit L2 norm)
%
% See bayes_warp_geom for the two representations.

T = size(Phi, 1);
dt = 1 / (T - 1);
V = Phi * c;
K = size(c, 2);
psi = zeros(T, K);
for k = 1:K
    v = V(:, k);
    if rep == 0
        r = sqrt(trapz(v.^2) * dt);
        if r < 1e-10
            p = ones(T, 1) + v;
        else
            p = cos(r) + sin(r) / r * v;
        end
        psi(:, k) = abs(p);
    else
        psi(:, k) = exp((v - max(v)) / 2);
    end
end
psi2 = psi.^2;
G = [zeros(1, K); cumsum(dt / 2 * (psi2(1:end-1, :) + psi2(2:end, :)), 1)];
GT = G(end, :);
gam = G ./ GT;
gam(end, :) = 1;
psi = psi ./ sqrt(GT);
end
