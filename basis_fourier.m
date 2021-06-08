function out = basis_fourier(f_domain, numBasis, fourier_p)
result = zeros(length(f_domain), 2*numBasis);
for i = 1:(2*numBasis)
    j = ceil(i/2);
    if (mod(i,2) == 1)
        result(:,i) = sqrt(2) * sin(2*j*pi*f_domain./fourier_p);
    end
    if (mod(i,2) == 0)
        result(:,i) = sqrt(2) * cos(2*j*pi*f_domain./fourier_p);
    end
    out.x = f_domain;
    out.matrix = result;
end