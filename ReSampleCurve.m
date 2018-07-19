function Xn = ReSampleCurve(X,N)
% RESAMPLECURVE Resample curve to have N points
% -------------------------------------------------------------------------
% Resample curve
% 
% Usage: Xn = ReSampleCurve(X,N)
%
% This function converts functions to srsf
%
% Input:
% X: matrix (nxT) of n dimensional curve with T sample points
% N: number of points
% 
% Output:
% X: matrix (nxN) of n dimensional curve with N sample points
[n,T] = size(X);
del = zeros(1,T);
del(1) = 0;
for r = 2:T
    del(r) = norm(X(:,r) - X(:,r-1));
end
cumdel = cumsum(del)/sum(del);

newdel = (0:N-1)/(N-1);

Xn = zeros(n,N);
for j=1:n
    Xn(j,:) = interp1(cumdel,X(j,1:T),newdel,'linear');
end