function p = q_to_curve(q,scale)
% Q_TO_CURVE Convert SRVF to curve
% -------------------------------------------------------------------------
% Convert to SRVF to curve
% 
% Usage: p = q_to_curve(q)
%
% This function converts SRVFs to curves
%
% Input:
% q: matrix (n,T) defining T points on n dimensional SRVF
% scale: scale of curve
% 
% Output:
% q: matrix of curve

[n,T] = size(q);
if nargin < 2
    scale = 1;
end

qnorm = zeros(1,T);
for i = 1:T
    qnorm(i) = norm(q(:,i),'fro');
end

p = zeros(n,T);
for i = 1:n
    p(i,:) = cumtrapz(q(i,:).*qnorm)/(T);
end

p = scale*p;
