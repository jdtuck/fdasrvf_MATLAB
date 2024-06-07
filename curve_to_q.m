function [q,len,lenq] = curve_to_q(p,closed)
% CURVE_TO_Q Convert curve to Square-Root Velocity Function
% -------------------------------------------------------------------------
% Convert to SRVF
% 
% Usage: q = curve_to_q(p)
%
% This function converts cruves to SRVFs
%
% Input:
% p: matrix (n,T) defining T points on n dimensional curve
% closed: if curve is closed
% 
% Output:
% q: matrix of SRVF
% len: length of curve
% lenq: length of SRVF

if nargin < 2
    closed = false;
end

[n,T] = size(p);
v=gradient(p,1/(T-1));


q=zeros(n,T);
len = sqrt(InnerProd_Q(sqrt(abs(v)),sqrt(abs(v))));
for i = 1:T
    L = sqrt(norm(v(:,i),'fro'));
    if L > 0.0001
        q(:,i) = v(:,i)/L;
    else
        q(:,i) = v(:,i)*0.0001;
    end
end

lenq = sqrt(InnerProd_Q(q,q));
q = q/lenq;

if closed
    q = ProjectC(q);
end
