function x_bound = mirrorbound_1(x, bidx, bdms)
% Calculate mirrored data off boundary of x, in dimension bdm
%   x(-k) = x(k) 
%   x(N+k) = x(N-k) for k > 0

    if not(iscell(bidx))
        bidx = {bidx};
    end

    dms = 1:ndims(x);

    % Construct multidimensional subscript
    inidx = cell(size(dms));

    for dm = dms
        inidx{dm} = ':';
    end

    M = size(x);
    
    for bi = 1:length(bdms)
        dm = bdms(bi);
        inidx{dm} = mirroridx(bidx{bi}, M(dm));
    end
    
    if(isvector(x) || isvector(inidx{1}))
        x_bound = x(inidx{:});
    else
        x_bound = x(sub2ind(size(x),inidx{:}));
    end