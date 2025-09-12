function H = Helmertsub(k)
% Helmertsub Helmert sub-matrix
%       Helmertsub(k) returns (k-1) x k Helmert sub-matrix
%
%

% Last updated Oct 20, 2009
% Sungkyu Jung

H = zeros(k-1,k);
for j= 1:(k-1)
    hj = -1/sqrt(j*(j+1));
    H(j,:) = [ones(1,j)*hj -j*hj zeros(1,k-j-1)];
end