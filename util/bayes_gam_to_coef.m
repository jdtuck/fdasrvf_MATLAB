function c = bayes_gam_to_coef(gam, Phi, rep)
% BAYES_GAM_TO_COEF Project warping functions onto warping basis coefficients
% -------------------------------------------------------------------------
% Inverse of bayes_coef_to_gam up to basis truncation.
%
% Usage: c = bayes_gam_to_coef(gam, Phi, rep)
%
% Input:
% gam: (T,K) warping functions on the uniform grid linspace(0,1,T)
% Phi: (T,D) orthonormal zero-mean basis (bayes_warp_basis)
% rep: 0 = shooting vector (v-transform), 1 = clr (h-transform)
%
% Output:
% c: (D,K) coefficients

[T, K] = size(gam);
dt = 1 / (T - 1);
w = dt * ones(T, 1);
w([1, T]) = dt / 2;
c = zeros(size(Phi, 2), K);
for k = 1:K
    g = (gam(:, k) - gam(1, k)) / (gam(end, k) - gam(1, k));
    gd = max(gradient(g, dt), 1e-8);
    if rep == 0
        psi = sqrt(gd);
        psi = psi / sqrt(sum(w .* psi.^2));
        v = inv_exp_map(ones(T, 1), psi);
        r = sqrt(sum(w .* v.^2));
        if r >= 0.95 * pi
            v = v * (0.95 * pi / r);
        end
    else
        v = log(gd);
        v = v - sum(w .* v);
    end
    c(:, k) = Phi.' * (w .* v);
end
end
