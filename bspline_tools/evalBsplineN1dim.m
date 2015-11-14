function f = evalBsplineN1dim( c, n, xs, bsh, boundaryfun )
% Evaluate B-spline across 1st non-singleton dimension
%   c - spline coefs
%   n - spline order
%   xs - points to evaluate at 1st non-singleton dimension
%   bsh - B-spline basis function shift
%   boundaryfun - boundary function used for padding the input

    if(nargin<5)
        boundaryfun = @mirrorbound_1;
    end

    m = round((n+1)/2);  
    
    % Remove leading singleton dimensions
    [c,lsdims] = shiftdim(c);
    
    fullidx = cell(size(1:ndims(c)));    
    
    for d = 1:ndims(c)
        fullidx{d} = ':';
    end
    
    % Pad columns
    c = boundaryfun(c, (-m+1):(size(c,1)+m), 1);
    
    csize = size(c);
    fsize = [length(xs), csize(2:end)];
    f = zeros(fsize);
    
    for xidx = 1:length(xs)
        
        x = xs(xidx);
        lmin = ceil(x+bsh-(n+1)/2);
        lmax = lmin + n;

        l = (lmin:lmax);

        cidx = fullidx;
        cidx{1} = l+m;
        
        outidx = fullidx;
        outidx{1} = xidx;
        
        f(outidx{:}) = sum(bsxfun(@times, c(cidx{:}), bsplineN(x - l + bsh, n)'), 1);
    end
    
    % Replace the leading singleton dimensions
    f = shiftdim(f, -lsdims);
end