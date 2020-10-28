function [dy, dx] = elastic_distance_curve(beta1, beta2, closed, method)
% ELASTIC_DISTANCE_CURVE Calculates the two elastic distances between two
% curves
% -------------------------------------------------------------------------
% This functions calculates the distances between curves,
% \eqn{D_y} and \eqn{D_x}, where curve 1 is aligned to curve 2
%
% Usage: [dy, dx] = elastic_distance_curve(beta1, beta2)
%        [dy, dx] = elastic_distance_curve(beta1, beta2, closed)
%
% Input:
% beta1: sample curve 1 (nxN matrix) n is dimension and N is number of
% amples
% beta2: sample curve 1
% closed: boolean if curve is closed (default = false)
% method: controls which optimization method (default="DP") options are
% Dynamic Programming ("DP") and Riemannian BFGS
% ("RBFGSM")
%
% Output
% dy: amplitude distance
% dx: phase distance
if nargin < 3
    closed = false;
    method = 'DP';
elseif nargin < 4
    method = 'DP';
end

N = size(beta1,2);
a=-calculateCentroid(beta1);
beta1 = beta1 + repmat(a,1,N);
a=-calculateCentroid(beta2);
beta2 = beta2 + repmat(a,1,N);
q1 = curve_to_q(beta1);
q2 = curve_to_q(beta2);

% Compute shooting vector from mu to q_i
[qn_t,~,gam] = Find_Rotation_and_Seed_unique(q1,q2,true,closed, method);
[qn_t,~] = Find_Best_Rotation(q1,qn_t);

q1dotq2=InnerProd_Q(q1,qn_t);

% Compute shooting vector
if q1dotq2>1
    q1dotq2=1;
end

dy = acos(q1dotq2);

time1 = linspace(0,1,N);
binsize = mean(diff(time1));
psi = sqrt(gradient(gam,binsize));
dx = acos(trapz(time1,psi));
