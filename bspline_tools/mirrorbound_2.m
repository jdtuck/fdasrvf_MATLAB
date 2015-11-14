function x_bound = mirrorbound_2(x, bidx, bdms)
    dms = 1:ndims(x);
    
    if not(iscell(bidx))
        bidx = {bidx};
    end

    % Construct multidimensional subscript
    inidx = cell(size(dms));

    for dm = dms
        inidx{dm} = ':';
    end

    M = size(x);
    
    for bi = 1:length(bdms)
        dm = bdms(bi);
        inidx{dm} = mirroridx_hs(bidx{bi}, M(dm));
    end
    
    x_bound = x(inidx{:});