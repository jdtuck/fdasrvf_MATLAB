function [q,len] = curve_to_q(p,scale,closed)
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
% scale: scale curve to length 1
% closed: if curve is closed
% 
% Output:
% q: matrix of SRVF
% len: length of curve

if nargin < 2
    scale = true;
end

if nargin < 3
    scale = true;
    closed = false;
end

[n,T] = size(p);
v=gradient(p,1/(T-1));


q=zeros(n,T);
for i = 1:T
    L = sqrt(norm(v(:,i),'fro'));
    if L > 0.0001
        q(:,i) = v(:,i)/L;
    else
        q(:,i) = v(:,i)*0.0001;
    end
end

len = sqrt(InnerProd_Q(q,q));
if scale
    q = q/len;
end

if closed
    q = ProjectC(q);
end