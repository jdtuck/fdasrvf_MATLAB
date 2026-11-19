function [nll, grad, SSE, gam, q2w] = bayes_warp_geom(c, Phi, rep, q1, q2, var1)
% BAYES_WARP_GEOM Negative log-likelihood and exact gradient of a warp
% -------------------------------------------------------------------------
% Model (SRVF space):  q1(t_k) = (q2 o gam)(t_k) sqrt(gam'(t_k)) + e_k,
% e_k ~ N(0, var1), with the warping function gam represented by basis
% coefficients c on the uniform grid t = linspace(0,1,T).
%
%   rep = 0 (shooting vector / v-transform): v = Phi*c is a tangent vector
%           at the identity of the Hilbert sphere, psi = cos|v| + sin|v|/|v| v,
%           gam = int psi^2. Requires |v| < pi (nll = Inf otherwise).
%   rep = 1 (centered log-ratio / h-transform): h = Phi*c,
%           gam' = exp(h) / int exp(h); unconstrained.
%
% Usage: [nll, grad, SSE, gam, q2w] = bayes_warp_geom(c, Phi, rep, q1, q2, var1)
%
% Input:
% c: (D,1) basis coefficients
% Phi: (T,D) basis matrix on the uniform grid (columns integrate to zero)
% rep: 0 = shooting vector, 1 = clr (h-transform)
% q1: (T,1) target SRVF
% q2: (T,1) SRVF to be warped
% var1: noise variance
%
% Batched use: c (D,N), q2 (T,N) and q1 (T,1) shared or (T,N) evaluate N
% independent warps; nll and SSE are (1,N), grad (D,N), gam and q2w (T,N).
%
% Output:
% nll: negative log-likelihood SSE/(2 var1) + T/2 log(var1)
% grad: (D,1) gradient of nll with respect to c
% SSE: sum of squared residuals
% gam: (T,1) warping function on the grid
% q2w: (T,1) warped SRVF (q2 o gam) sqrt(gam')
%
% A compiled version (bayes_warp_geom_c) is used automatically when built
% (buildfile task mex:bayes_warp_geom).

persistent has_mex
if isempty(has_mex)
    has_mex = exist('bayes_warp_geom_c', 'file') == 3;
end
if has_mex
    [nll, grad, SSE, gam, q2w] = bayes_warp_geom_c(c, Phi, double(rep), q1, q2, var1);
    return
end

if size(c, 2) > 1
    % batched: one independent warp per column
    N = size(c, 2);
    T = size(Phi, 1);
    nll = zeros(1, N); grad = zeros(size(c)); SSE = zeros(1, N);
    gam = zeros(T, N); q2w = zeros(T, N);
    for i = 1:N
        q1i = q1(:, min(i, size(q1, 2)));
        [nll(i), grad(:,i), SSE(i), gam(:,i), q2w(:,i)] = ...
            bayes_warp_geom(c(:,i), Phi, rep, q1i, q2(:,i), var1);
    end
    return
end

T = size(Phi, 1);
D = size(Phi, 2);
dt = 1 / (T - 1);
w = dt * ones(T, 1);
w([1, T]) = dt / 2;
v = Phi * c;

if rep == 0
    r = sqrt(sum(w .* v.^2));
    if r >= pi
        nll = Inf; grad = zeros(D, 1); SSE = Inf;
        gam = linspace(0, 1, T).'; q2w = q2;
        return
    end
    if r < 1e-10
        s = 1; ds = 0; sr = 0;
        psi = 1 + v;
    else
        s = sin(r) / r;
        ds = (r * cos(r) - sin(r)) / r^2;
        sr = sin(r);
        psi = cos(r) + s * v;
    end
    sgn = sign(psi);
    sgn(sgn == 0) = 1;
    apsi = abs(psi);
else
    psi = exp((v - max(v)) / 2);
    sgn = ones(T, 1);
    apsi = psi;
end

% gam = G / G(T), G = cumtrapz(psi^2)
psi2 = psi.^2;
G = [0; cumsum(dt / 2 * (psi2(1:end-1) + psi2(2:end)))];
GT = G(end);
gam = G / GT;
gam(end) = 1;

% linear interpolation of q2 at gam on the uniform grid; q2dg is the slope
% of the interpolant, so the gradient is exact for the discretized model
[q2g, q2dg] = interp_uniform(q2, gam, T);

a = 1 / sqrt(GT);
q2w = q2g .* apsi * a;
res = q1 - q2w;
SSE = sum(res.^2);
nll = SSE / (2 * var1) + T / 2 * log(var1);

if nargout > 1
    % d(sum res.*q2w)/dc written as beta' * dpsi/dc: the cumulative
    % trapezoid in gam is reverse-accumulated so no (T,D) arrays are formed
    wk = res .* a;
    A = wk .* q2dg .* apsi;          % weight on dgam
    Bv = wk .* q2g .* sgn;           % weight on dpsi
    kappa = sum(res .* q2w) / (2 * GT);
    S = flipud(cumsum(flipud(A)));   % suffix sums of A
    coefP = zeros(T, 1);
    coefP(1:T-1) = dt / 2 * S(2:T);
    coefP(2:T) = coefP(2:T) + dt / 2 * S(2:T);
    alpha = coefP / GT - (sum(A .* gam) / GT + kappa) * w; % weight on d(psi^2)
    beta = 2 * alpha .* psi + Bv;
    if rep == 0
        gvec = s * (Phi.' * beta);
        if r >= 1e-10
            dr = (Phi.' * (w .* v)) / r;
            gvec = gvec + (-sr * sum(beta) + ds * sum(beta .* v)) * dr;
        end
    else
        gvec = 0.5 * (Phi.' * (beta .* psi));
    end
    grad = -gvec / var1;
end
end

function [yi, ydi] = interp_uniform(y, x, T)
pos = x * (T - 1) + 1;
pos = min(max(pos, 1), T);
k = min(floor(pos), T - 1);
fr = pos - k;
yi = y(k) .* (1 - fr) + y(k + 1) .* fr;
ydi = (y(k + 1) - y(k)) * (T - 1);
end
