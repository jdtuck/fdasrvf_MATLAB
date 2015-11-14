function y = filterAPrev1ord( x, zi )
% filterAPrev1ord - reverse direction 1st order allpole filter
%
%   cp = filterAPrev1ord( x, zi ) filters input signal x 
%   according to the transfer function
%   H(z) = -zi / (1 - zi * z)
%
%   The reflected input data is used as terminal conditions
%   x(N+k) = x(N-k) for k > 0
%
%   The filter operates along the first non-singleton dimension of x
%
%   This function requires the Signal Processing toolbox
%
%   See also filterAPfwd1ord, filterAPsym2ord

b1 = -zi;
a1 = [1, -zi];


dms = 1:ndims(x);
% First non-singleton dimension
fnsdm = find(size(x)>1, 1);

% We will filter in the forward direction, flip
x = flipdim(x, fnsdm);
y = zeros(size(x));


% Construct input and output subscripts
fullidx = cell(size(dms));
for dm = dms;
    fullidx{dm} = ':';
end

xidx1 = fullidx;
xidx1{fnsdm} = 1;
xidx2 = fullidx;
xidx2{fnsdm} = 2;

yidx1 = fullidx;
yidx1{fnsdm} = 1;

% Final sample calculated using a closed-form expression
y(yidx1{:}) = (zi/(zi*zi-1)) * (x(xidx1{:})+zi*x(xidx2{:}));

% Calculate state variables for the filter function
filtic_b1_a1 = @(y,x) filtic(b1,a1,y,x);
zn0 = arrayfun(filtic_b1_a1, y(yidx1{:}), x(xidx1{:}));

% Remaining samples calculated using the filter function
xidxt = fullidx;
xidxt{fnsdm} = 2:size(x, fnsdm);
y(xidxt{:}) = filter(b1, a1, x(xidxt{:}), zn0);
y = flipdim(y, fnsdm);

end

