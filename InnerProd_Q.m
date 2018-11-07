function val = InnerProd_Q(q1,q2)
% INNERPROD_Q Calculate inner product of two SRVFs
% -------------------------------------------------------------------------
% Convert to SRSF
% 
% Usage: val = InnerProd_Q(q1,q2)
%
% This function calculates the innter product between two curves
%
% Input:
% q1: matrix of (nxT) of SRVF1
% q2: matrix of (nxT) of SRVF2
% 
% Output:
% val: inner product
[~,T] = size(q1);
val=sum(sum(q1.*q2))/T;

return;