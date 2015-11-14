function c = filterAPsym2ord( x, zi, boundaryfun )
% filterAPsym2ord - symmetrical 2. order allpole filter
%   c = filterAPsym2ord( x, zi) filters signal x according to
%   the transfer function 
%   H(z) = - zi / (1 - zi * z)*(1 - zi * z^-1) 
%
%   The boundary signal is calculated by the function given as
%   the third argument:
%   c = filterAPsym2ord( x, zi, boundaryfun)
%
%   The default boundary function is mirrorbound_1
%
%   The filter operates along the first non-singleton dimension of x
%
%   This function requires the Signal Processing toolbox
%
%   See also mirrorbound_1

b0 = 1;
a0 = [1, -zi];

b1 = -zi;
a1 = [1, -zi];

if(nargin<3)
    boundaryfun = @mirrorbound_1;
end

fnsdim = find(size(x)>1, 1); % first non-singleton dimension

% Number of reflected data points used to
% initialize the IIR filter state
k0 = ceil(log(eps) / log(abs(zi)));
%k0 = min(k0,size(x,fnsdim)-1);

x_bo = boundaryfun(x, -(k0-1):0, fnsdim);
[cp_thrash, zp1] = filter(b0, a0, x_bo);

[cp, zp2] = filter(b0, a0, x, zp1);

% Calculate boundary on the upper end
x_bo_p = boundaryfun(x, (1:k0) + size(x,fnsdim), fnsdim);
cp_bo = filter(b0, a0, x_bo_p, zp2);

% We will filter in the reverse direction, flip
cp = flipdim(cp, fnsdim);
cp_bo = flipdim(cp_bo, fnsdim);

[cp_thrash, zn0] = filter(b1, a1, cp_bo);

c = filter(b1, a1, cp, zn0);
c = flipdim(c, fnsdim);


end

