function f = evalBsplineNpoint( c, n, p, bsh )
% Evaluate centered spline at single point
%   c - spline coefs
%   n - spline order
%   m - padding
%   p - point to evaluate
%   bsh - basis shift
   
    dms = length(p);
    bspl = cell(size(p));
    cidx = cell(size(p));

    for dm = 1:dms
        x = p(dm);
        lmin = ceil(x+bsh(dm)-(n(dm)+1)/2);
        lmax = lmin + n(dm);

        l = (lmin:lmax);

        cidx{dm} = mirroridx(l, size(c,dm));

        bspl{dm} = shiftdim( bsplineN(x -l + bsh(dm), n(dm)), -dm + 2);
    end

    cblk = c(cidx{:});

    for dm = 1:dms
        % Since the B-spline basis is separable, we multiply with it
        % for each dimension
        cblk = bsxfun(@times, cblk, bspl{dm});
    end

    f = sum(cblk(:));
end