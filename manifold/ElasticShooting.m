function q2n = ElasticShooting(q1,v)
% ElasticShooting Calculate Shooting Vector from v to q1
% -------------------------------------------------------------------------
% 
%
% Usage: q2n = ElasticShooting(q1,v)
%
% This function computes the shooting of q1 to v
%
% Input:
% q1: vector of srvf
% v: shooting vector
%
% Output:
% q2n: srvf

d = sqrt(InnerProd_Q(v,v));
if d < 0.00001
    q2n = q1;
else
    q2n = cos(d)*q1 + (sin(d)/d)*v;
    q2n = ProjectC(q2n);
end