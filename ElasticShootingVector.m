function [v,d,q2n] = ElasticShootingVector(q1,q2,reparamFlag)
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
%
% Output:
% v: shooting vector
% d: distance
% q2n: aligned srvf


q2n = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag,1);
q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));

d = acos(InnerProd_Q(q1,q2n));
if d < 0.0001
    v = zeros(size(q1));
else
    v = (d/sin(d))*(q2n - cos(d)*q1);
end