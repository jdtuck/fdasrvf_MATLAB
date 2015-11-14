A = imread('lena512.bmp');
A_sz = size(A);
A_max = double(max(A(:)));
A_min = double(min(A(:)));

A = A(1:(floor(A_sz(1)/2)*2),1:(floor(A_sz(2)/2)*2));

Ab = Bspline(double(A), 3);

Abr = reduce(Ab);
Abe = expand(Abr);

Abr_double = double(Abr);
Ab_double = double(Ab);
Abe_double = double(Abe);
10*log10(numel(A)*(A_max-A_min).^2 ./ sum(sum((Abe_double-Ab_double).^2)))

imshow(Abe_double/256);
