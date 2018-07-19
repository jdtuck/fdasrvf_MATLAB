function val = InnerProd_Q(q1,q2)
% INNERPROD_Q Calculate inner product of two SRVFs
% -------------------------------------------------------------------------
% Convert to SRSF
% 
% Usage: val = InnerProd_Q(q1,q2)
%
% This function converts functions to srsf
%
% Input:
% q1: matrix of (nxT) of SRVF1
% q2: matrix of (nxT) of SRVF2
% 
% Output:
% val: inner product
[~,T] = size(q1);
val = trapz(linspace(0,1,T),sum(q1.*q2));

return;