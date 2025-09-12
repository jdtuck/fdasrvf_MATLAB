function Logpx = LogNPd(x)
% LOGNP Riemannian log map at North pole of S^k
%       LogNP(x) returns k x n matrix where each column is a point on tangent
%       space at north pole and the input x is (k+1) x n matrix where each column 
%       is a point on a sphere.
%
%
%   See also ExpNPd.

% Last updated Oct 20, 2009
% Sungkyu Jung

[d, ~] = size(x);
scale = acos(x(end,:))./sqrt(1-x(end,:).^2);
scale(isnan(scale)) =1;
Logpx = repmat(scale,d-1,1).*x(1:(end-1),:);


