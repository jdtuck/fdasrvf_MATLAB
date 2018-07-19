function q = curve_to_q(p)
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
% 
% Output:
% q: matrix of SRVF
[n,N] = size(p);
v = zeros(n,N);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N));
end

len = sum(sqrt(sum(v.*v)))/N;
v = v/len;
L = zeros(1,N);
q = zeros(n,N);
for i = 1:N
    L(i) = sqrt(norm(v(:,i)));
    if L(i) > 0.0001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = v(:,i)*0.0001;
    end
end