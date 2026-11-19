function [c, nll, g, aux, accept, nsteps] = infhmc_step(c, nll, g, aux, geom, C, epsilon, L)
% INFHMC_STEP One preconditioned (infinite-dimensional) HMC transition
% -------------------------------------------------------------------------
% Samples pi(c) ~ exp(-nll(c)) N(c; 0, diag(C)) with the splitting
% integrator of Beskos et al. (2011): the Gaussian prior is integrated
% exactly by a rotation and the likelihood enters through half kicks
% v <- v - epsilon/2 C grad(nll). The integrator is reversible and volume
% preserving, so acceptance uses the exact change in
% H = nll(c) + c'C^{-1}c/2 + v'C^{-1}v/2.
%
% Usage: [c, nll, g, aux, accept] = infhmc_step(c, nll, g, aux, geom, C, epsilon, L)
%
% Input:
% c: (D,1) current state
% nll, g, aux: negative log-likelihood, its gradient and any extra output
%              of geom at c (must be evaluated with the current parameters)
% geom: function handle [nll, g, aux] = geom(c)
% C: (D,1) prior variances
% epsilon: step size (rotation angle)
% L: maximum number of leapfrog steps (number used is uniform on 1..L)
%
% Output:
% c, nll, g, aux: new state and its geometry
% accept: true if the proposal was accepted
% nsteps: number of leapfrog steps used

C = C(:);
v = sqrt(C) .* randn(size(c));
H0 = nll + 0.5 * sum(c.^2 ./ C) + 0.5 * sum(v.^2 ./ C);

nsteps = ceil(rand * L);
q = c;
gq = g;
cs = cos(epsilon);
sn = sin(epsilon);
nll_q = nll;
aux_q = aux;
for l = 1:nsteps
    v = v - epsilon / 2 * C .* gq;
    q_new = q * cs + v * sn;
    v = v * cs - q * sn;
    q = q_new;
    [nll_q, gq, aux_q] = geom(q);
    if ~isfinite(nll_q)
        break
    end
    v = v - epsilon / 2 * C .* gq;
end

H1 = nll_q + 0.5 * sum(q.^2 ./ C) + 0.5 * sum(v.^2 ./ C);
logRatio = H0 - H1;
accept = isfinite(logRatio) && (log(rand) < logRatio);
if accept
    c = q;
    nll = nll_q;
    g = gq;
    aux = aux_q;
end
end
