function [u2N, c] = u2N_FIR_coefs(N)
% u2n_FIR_coefs - filter coefficients for B-spline scale relation
%
% [u2N, c] = u2N_FIR_coefs(N) generates the FIR filter u2N which relates two
% scalings of the b-spline basis:
% 
% b_N_2 = c*conv(u2N, b_N)
%
% where b_N and b_N_2 are two discretized versions of the basic B-spline
% of order N, b_N_2 is dilated with a factor of two:
%
%       b_N_2 = beta_N(k/2) for k = -Inf, Inf
%       b_N = beta_N(k) for k = -Inf, Inf
%
% c is an additional scale factor.


u2N = zeros(1,N+2);
c = 2^(-N);

for k = (-(N+1)/2):((N+1)/2);
    u2N(k+(N+1)/2+1) = nchoosek((N+1),k+(N+1)/2);
end
    