function [v,d,q2n] = ElasticShootingVector(q1,q2,reparamFlag,scale)
% ElasticShootingVector Find Shooting Vector between to SRVFs
% -------------------------------------------------------------------------
% Compute Shooting for two SRVFs
%
% Usage: [v,d,q2n] = ElasticShootingVector(q1,q2,reparamFlag)
%
% This function computes the shooting between two SRVFS
%
% Input:
% q1: srvf of curve 1
% q2: srvf of curve 2
% reparamFlag: compute reparameterization
% scale: scale of q2 (default: 0. or don't include)
%
% Output:
% v: shooting vector
% d: distance
% q2n: aligned srvf

if nargin < 4
    scale = 0;
end
if scale == 0
    ifscale = false;
end

q2n = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag,1);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));

d = acos(InnerProd_Q(q1,q2n));
if d < 0.0001
    v = zeros(size(q1));
else
    v = (d/sin(d))*(q2n - cos(d)*q1);
end

if ifscale
    v = [v(:); scale];
end
