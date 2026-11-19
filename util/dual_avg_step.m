function [epsilon, state] = dual_avg_step(state, alpha, target)
% DUAL_AVG_STEP Dual-averaging step-size adaptation for HMC
% -------------------------------------------------------------------------
% Nesterov dual averaging as in Hoffman and Gelman (2014), "The No-U-Turn
% Sampler", Algorithm 5. Call once per burn-in iteration with the acceptance
% probability of the transition; at the end of burn-in use
% exp(state.log_eps_bar), the averaged step size, which is much less noisy
% than the last iterate. Works elementwise for a vector of step sizes.
%
% Usage: state = dual_avg_step(epsilon0)            % initialize
%        [epsilon, state] = dual_avg_step(state, alpha, target)
%
% Input:
% epsilon0: initial step size(s)
% state: adaptation state
% alpha: acceptance probability(ies) of the last transition, in [0,1]
% target: target acceptance probability
%
% Output:
% epsilon: step size(s) for the next iteration (capped at pi/2)
% state: updated state; state.log_eps_bar holds the averaged log step size

if nargin == 1
    eps0 = state;
    epsilon = struct('mu', log(10 * eps0), 'Hbar', zeros(size(eps0)), ...
        'log_eps_bar', log(eps0), 'n', 0);
    return
end

t0 = 10;
gamma = 0.05;
kappa = 0.75;

state.n = state.n + 1;
n = state.n;
alpha = reshape(alpha, size(state.Hbar));
state.Hbar = (1 - 1 / (n + t0)) * state.Hbar + (target - alpha) / (n + t0);
log_eps = min(state.mu - sqrt(n) / gamma * state.Hbar, log(pi / 2));
eta = n^(-kappa);
state.log_eps_bar = eta * log_eps + (1 - eta) * state.log_eps_bar;
epsilon = exp(log_eps);
end
