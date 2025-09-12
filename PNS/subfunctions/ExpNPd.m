function Exppx = ExpNPd(x)
% EXPNP Riemannian exponential map at North pole of S^k
%       ExpNP(v) returns (k+1) x n matrix where each column is a point on a 
%                sphere and the input v is k x n matrix where each column 
%                is a point on tangent  space at north pole.
%
%
%   See also LogNPd.

% Last updated Oct 20, 2009
% Sungkyu Jung

[d, ~] = size(x);
nv = sqrt(sum(x.^2));
Exppx = [repmat(sin(nv)./nv,d,1).*x ;cos(nv)];
Exppx(:,nv < 1e-16) = repmat([zeros(d,1);1],1,sum(nv < 1e-16));
