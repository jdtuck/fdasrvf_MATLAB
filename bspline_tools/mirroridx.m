function idx = mirroridx(idx,M)
    % Mirror an index value around 1 and M
    
    idx = abs1(mod(idx + M - 3, 2*M - 2) - M + 3);
end

function y = abs1(x)
    y = abs(x-1) + 1;
end