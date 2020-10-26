function Xn = ReSampleCurve(X,N,closed)
% RESAMPLECURVE Resample curve to have N points
% -------------------------------------------------------------------------
% Resample curve
%
% Usage: Xn = ReSampleCurve(X,N)
%
% This function resamples a curve on N points
%
% Input:
% X: matrix (nxT) of n dimensional curve with T sample points
% N: number of points
% closed: if curve is closed (default T)
%
% Output:
% X: matrix (nxN) of n dimensional curve with N sample points
if nargin < 3
    closed = false;
end
[n,T] = size(X);
del = zeros(1,T);
del(1) = 0;
tst = X(:,2)-X(:,1);
if tst(1) < 0
    X = fliplr(X);
end
for r = 2:T
    del(r) = norm(X(:,r) - X(:,r-1));
end
cumdel = cumsum(del)/sum(del);

newdel = linspace(0,1,N);

Xn = zeros(n,N);
for j=1:n
    Xn(j,:) = interp1(cumdel,X(j,:),newdel,'makima');
end

if closed
    q = curve_to_q(Xn);
    qn = ProjectC(q);
    Xn = q_to_curve(qn);
end
