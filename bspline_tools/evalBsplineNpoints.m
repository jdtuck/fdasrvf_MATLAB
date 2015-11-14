function f = evalBsplineNpoints( c, n, p, bsh, boundaryfun )
% Evaluate polynomial B-spline at an array of points
%   c - spline coefs
%   n - spline order
%   p - points to evaluate
%   bsh - basis shift
%   boundaryfun - boundary function for signal extension
    if(nargin<5)
        boundaryfun = @mirrorbound_1;
    end   


    dms = length(p);

    % Construct a grid of the neighborhood
    % to evaluate around each point
    nvecs = cell(1,dms);
    for dm = 1:dms
        nvecs{dm} = 0:n(dm);
    end
    lgrid = cell(1,dms);
    [lgrid{:}] = ndgrid(nvecs{:});
    
    csum = zeros(size(p{1}));
    
    % Loop through neighborhood
    for lidx = 1:prod((n+1))
        bspl = ones(size(p{1}));
        cidx = cell(size(p));
        % Calculate index and Bspline function for each dimension
        for dm = 1:dms
            x = p{dm};
            lmin = ceil(x+bsh(dm)-(n(dm)+1)/2);
            l = lmin + lgrid{dm}(lidx);
            cidx{dm} = l;

            bspl = bspl.*bsplineN(x -l + bsh(dm), n(dm));
        end

        % Accumulate B-spline in output domain
        csum = csum + boundaryfun(c, cidx, 1:dms).*bspl;
    end

    f = csum;
end