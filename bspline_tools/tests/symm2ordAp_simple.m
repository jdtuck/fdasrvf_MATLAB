function c = symm2ordAp_simple( x, zi )
% symmetrical 2. order allpole filter
%   x  - input signal
%   zi - pole of transfer function


cp = forward1ordAp_simple(x,zi);

c = reverse1ordAp_simple(cp, zi);


end

