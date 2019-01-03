function basismat = create_basis_legendre(m, N)
% CREATE_BASIS_LEGENDRE Create Basis Matrix for Legendre polynomials
%%% compute the Legendre polynomial coefficients matrix coeff
% coeff(i,j) gives the polynomial coefficient for term x^{j-1} in P_{i-1}(x)
if N > 1
    coeff = zeros(N+1);
    coeff([1 N+3]) = 1; % set coefficients of P_0(x) and P_1(x)
    % now compute for higher order: nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    for ii = 3:N+1
        coeff(ii,:) = (2-1/(ii-1))*coeff(ii-1,[end 1:end-1]) - (1-1/(ii-1))*coeff(ii-2,:);
    end
else
    % simple case
    coeff = eye(N+1);
end
x = linspace(-1,1,m);%-1:2/(m-1):1; % Legendre polynomials are supported for |x|<=1
x = x(:);
%%% Evaluate the polynomials for every element in X
basismat = cumprod([ones(m,1) x(:,ones(1,N))], 2) * coeff.';

% A = (D.'*D)\(D.'*Y);
% Y2 = D*A; 