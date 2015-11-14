function cp = filterAPfwd1ord( x, zi )
% filterAPfwd1ord - forward direction 1st order allpole filter
%
%   cp = filterAPfwd1ord( x, zi ) filters input signal x
%   according to the transfer function
%   H(z) = 1 / (1 - zi * z^-1)
%
%   The reflected input data is used as initial conditions:
%   x(-k) = x(k) 
%
%   The filter operates along the first non-singleton dimension of x
%
%   This function requires the Signal Processing toolbox
%   See also filterAPrev1ord, filterAPsym2ord

b0 = 1;
a0 = [1, -zi];

cp = zeros(size(x));

dms = 1:ndims(x);
fnsdm = find(size(x)>1, 1); % first non-singleton dimension

% Construct multidimensional subscripts
inidx = cell(size(dms));
oidx = cell(size(dms));

for dm = dms
    inidx{dm} = ':';
    oidx{dm} = ':';
end

% Number of reflected data points used to calculate the first sample
% and initialize the IIR filter state
k0 = round(log(0.5*sqrt(eps)) / log(abs(zi)));
k0 = min(k0,size(x,fnsdm)-1);


% Calculate the first output sample,
% by excplicitly calculating the first k0+1 IIR fiter coefs
% and applying to reflected input data (all along fns dimension)
inidx{fnsdm} = 1:(k0+1);
oidx{fnsdm} = 1;

zpows = cumprod(zi*ones(1,k0))';
zpows = [1; zpows];
zpows = shiftdim(zpows, -(fnsdm-1));
cp(oidx{:}) = sum(bsxfun(@times, zpows, x(inidx{:})), fnsdm);

% Calculate state values to use in the filter function 
inidx{fnsdm} = 1;
filtic_b0_a0 = @(y,x) filtic(b0,a0,y,x);
zp1 = arrayfun(filtic_b0_a0, cp(inidx{:}), x(inidx{:}) );

% Calculate the remaining output samples using SPTOOLBOX filter function
oidx{fnsdm} = 2:size(cp, fnsdm);
inidx{fnsdm} = 2:size(x, fnsdm);
cp(oidx{:}) = filter(b0, a0, x(inidx{:}), zp1);

end

