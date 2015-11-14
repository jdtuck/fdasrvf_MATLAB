function smooth_data(f, sparam)
[M,N] = size(f);
for r = 1:sparam
    for i = 1:N
        f(2:(M-1),i) = (f(1:(M-2),i)+2*f(2:(M-1),i) + f(3:M,i))/4;
    end
end