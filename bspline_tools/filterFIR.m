function cp = filterFIR( b, x, boundaryfun )
% filterFIR - FIR filter with functional boundary conditions
%
%   The filter operates along the first non-singleton dimension of x
%
%   This function requires the Signal Processing toolbox
%   See also filterAPrev1ord, filterAPsym2ord

if(nargin<3)
    boundaryfun = @mirrorbound_1;
end

dms = 1:ndims(x);
fnsdm = find(size(x)>1, 1); % first non-singleton dimension


% Number of data points used to extend the input
m = (length(b)-1)/2;

% Extend the input using the boundary function
z = boundaryfun(x, (-m+1):(size(x,fnsdm)+m), fnsdm);

% Calculate the output samples using SPTOOLBOX filter function
cp = filter(b, 1, z);

% Construct multidimensional subscripts to remove pad from output
oidx = cell(size(dms));

for dm = dms
    oidx{dm} = ':';
end 
oidx{fnsdm} = (m+m+1):(size(z, fnsdm));

cp = cp(oidx{:});

end

