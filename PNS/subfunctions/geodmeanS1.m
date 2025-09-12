function [geodmean, geodvar]= geodmeanS1(theta)
% geodmeanS1    geodesic mean of data on S^1 (Circle) by S. Lu and V. Kulkarni
%               method - gives all multiples of geodesic mean set.
%
% Inputs:   theta   - a column vector of angles
%                 
% Output:
%           geodmean    - geodesic mean on S^1
%           geodvar     - geodesic variance on S^2
% 
%    See also geodmeanS1mat, geoddistS2, geodmeanS2, geodmeanS2grid, geodmeanS2q, geodmeanS1.

% 2008.07.08, last updated Aug 14, 2009
% Sungkyu Jung



n = length(theta);
meancandi = mod(mean(theta)+2*pi*(0:n-1)/n,2*pi);
theta = mod(theta,2*pi);
geodvar = zeros(n,1);
for i=1:n
    v= meancandi(i);
    dist2 = min([(theta-v).^2,(theta-v+2*pi).^2,  (v-theta+2*pi).^2],[],2);
    geodvar(i) = sum(dist2);
end
[~, ind] = min(geodvar);
geodmean = mod(meancandi(ind),2*pi);
geodvar = geodvar(ind)/n;
