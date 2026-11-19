function [c, nll, g, aux, accept, alpha] = infhmc_step_batch(c, nll, g, aux, geom, C, epsilon, L)
% INFHMC_STEP_BATCH Preconditioned HMC transitions for independent blocks
% -------------------------------------------------------------------------
% Runs infhmc_step on N conditionally independent blocks at once (one
% column of c per block), sharing a target of the form
% prod_i exp(-nll_i(c_i)) N(c_i; 0, diag(C)). All blocks integrate the same
% (random) number of leapfrog steps, with their own step sizes, and each is
% accepted or rejected on its own energy change; since the blocks do not
% interact, every column is a valid HMC transition for its own conditional
% posterior. geom is evaluated once per leapfrog step for all blocks, so a
% batched likelihood (e.g. bayes_warp_geom with matrix input) amortizes
% call overhead.
%
% Usage: [c, nll, g, aux, accept] = infhmc_step_batch(c, nll, g, aux, geom, C, epsilon, L)
%
% Input:
% c: (D,N) current states
% nll: (1,N) negative log-likelihoods at c
% g: (D,N) gradients of nll at c
% aux: (P,N) extra output of geom at c (may be empty)
% geom: function handle [nll, g, aux] = geom(c) for all N columns
% C: (D,1) prior variances
% epsilon: (N,1) or scalar step sizes
% L: maximum number of leapfrog steps (number used is uniform on 1..L)
%
% Output:
% c, nll, g, aux: new states and their geometry
% accept: (N,1) logical, true where the proposal was accepted
% alpha: (N,1) acceptance probabilities min(1, exp(-dH)) (0 where the
%        trajectory left the support), for step-size adaptation

C = C(:);
[D, N] = size(c);
epsilon = reshape(epsilon, 1, []) .* ones(1, N);
v = sqrt(C) .* randn(D, N);
H0 = nll + 0.5 * sum(c.^2 ./ C, 1) + 0.5 * sum(v.^2 ./ C, 1);

nsteps = ceil(rand * L);
q = c;
gq = g;
cs = cos(epsilon);
sn = sin(epsilon);
nll_q = nll;
aux_q = aux;
failed = false(1, N);
for l = 1:nsteps
    v = v - (epsilon / 2) .* C .* gq;
    q_new = q .* cs + v .* sn;
    v = v .* cs - q .* sn;
    q = q_new;
    [nll_q, gq, aux_q] = geom(q);
    % a block whose trajectory leaves the support (nll = Inf) is rejected,
    % as in infhmc_step; zero its force so it stays finite
    bad = ~isfinite(nll_q);
    failed = failed | bad;
    gq(:, bad) = 0;
    v = v - (epsilon / 2) .* C .* gq;
end

H1 = nll_q + 0.5 * sum(q.^2 ./ C, 1) + 0.5 * sum(v.^2 ./ C, 1);
logRatio = H0 - H1;
ok = ~failed & isfinite(logRatio);
alpha = zeros(N, 1);
alpha(ok) = min(1, exp(logRatio(ok)));
accept = (ok & (log(rand(1, N)) < logRatio)).';
c(:, accept) = q(:, accept);
nll(accept) = nll_q(accept);
g(:, accept) = gq(:, accept);
if ~isempty(aux)
    aux(:, accept) = aux_q(:, accept);
end
end
