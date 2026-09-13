function out = basis_fourierd(f_domain, numBasis)
% BASIS_FOURIERD Fourier basis evaluated on a grid
% -------------------------------------------------------------------------
% Builds the 2*NUMBASIS Fourier basis functions
% sin(2*j*pi*x)/sqrt(pi) and cos(2*j*pi*x)/sqrt(pi), j = 1, ..., NUMBASIS,
% interleaved as [sin_1, cos_1, sin_2, cos_2, ...].
%
% Usage: out = basis_fourierd(f_domain, numBasis)
%
% Input:
% f_domain: vector of grid points
% numBasis: number of sine/cosine pairs
%
% Output:
% out.x: the grid, as a column vector
% out.matrix: (length(f_domain) x 2*numBasis) matrix of basis values

f_domain = f_domain(:);
arg = 2 * pi * f_domain * (1:numBasis);

result = zeros(numel(f_domain), 2 * numBasis);
result(:, 1:2:end) = sin(arg) / sqrt(pi);
result(:, 2:2:end) = cos(arg) / sqrt(pi);

out.x = f_domain;
out.matrix = result;
end
