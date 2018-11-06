function [q,len] = curve_to_q(p,scale)
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
% 
% Output:
% q: matrix of SRVF
% len: length of curve

if nargin < 2
    scale = true;
end

[n,T] = size(p);
v=gradient(p,1/(T-1));

len = sum(sqrt(sum(v.*v)))/T;
if scale
    v = v/len;
end
q=zeros(n,T);
for i = 1:T
    L = sqrt(norm(v(:,i)));
    if L > 0.0001
        q(:,i) = v(:,i)/L;
    else
        q(:,i) = v(:,i)*0.0001;
    end
end