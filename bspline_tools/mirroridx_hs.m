function idx0 = mirroridx_hs(idx,M)
    % Mirror an index value around 0.5 and M+0.5
    
    idx = idx-1;
    
    idx0 = zeros(size(idx));
    idx_pos = idx>=0;
    idx_neg = idx<0;
    idx0(idx_pos) = mod(idx(idx_pos), 2*M);
    idx0(idx_neg) = mod(2*M-1-idx(idx_neg), 2*M);
    
    idx0(idx0>=M) = -idx0(idx0>=M)-1+2*M;
    
    idx0 = idx0+1;
end