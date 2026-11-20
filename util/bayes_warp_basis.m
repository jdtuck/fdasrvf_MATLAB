function Phi = bayes_warp_basis(t, ncoef, type)
% BAYES_WARP_BASIS Zero-mean orthonormal basis for warping coefficients
% -------------------------------------------------------------------------
% Returns ncoef basis functions on t in [0,1] that integrate to zero (so
% they lie in the tangent space at the identity / give a centered h) and
% are orthonormal in L2[0,1].
%
% Usage: Phi = bayes_warp_basis(t, ncoef, type)
%
% Input:
% t: grid on [0,1]
% ncoef: number of basis functions (even for 'fourier')
% type: 'fourier' (sqrt(2) sin(2 pi j t), sqrt(2) cos(2 pi j t)) or
%       'legendre' (shifted Legendre polynomials of degree 1..ncoef)
%
% Output:
% Phi: (length(t), ncoef) basis matrix

t = t(:);
Phi = zeros(length(t), ncoef);
if strcmpi(type, 'fourier')
    for i = 1:ncoef
        j = ceil(i/2);
        if mod(i, 2) == 1
            Phi(:, i) = sqrt(2) * sin(2*j*pi*t);
        else
            Phi(:, i) = sqrt(2) * cos(2*j*pi*t);
        end
    end
elseif strcmpi(type, 'legendre')
    x = 2*t - 1;
    for i = 1:ncoef
        tmp = legendre(i, x);
        Phi(:, i) = sqrt(2*i + 1) * tmp(1, :).';
    end
else
    error('Incorrect Basis Specified')
end
end
