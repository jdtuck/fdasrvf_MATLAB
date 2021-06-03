function out = basis_fourierd(f_domain, numBasis)
% construct fourier basis
result = zeros(length(f_domain), 2*numBasis);
for i = 1:(2*numBasis)
    j = ceil(i/2);
    if (mod(i,2) == 1)
        result(:,i) = 1/sqrt(pi) * sin(2*j*pi*f_domain);
    end
    if (mod(i,2) == 0)
        result(:,i) = 1/sqrt(pi) * cos(2*j*pi*f_domain);
    end
    out.x = f_domain;
    out.matrix = result;
end