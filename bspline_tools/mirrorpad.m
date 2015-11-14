function y = mirrorpad(x, m)
    % Pad the array x with a mirror image of itself
    % If m is a scalar, pads the first non-singleton dimension of x with m
    % elements
    % When m is a vector, pads each dimension according to elements of m.

    insize = size(x);
    
    dms = ndims(x);
    
    % Build a subscripting cell array programmatically
    
    idcell = cell(size(insize));
    
    if(length(m)==1)
        fns = 0; % Found non-singleton dimension
        for dm = 1:dms
            if(insize(dm)>1 && fns==0)
                idcell{dm} = mirroridx((-m+1):(insize(dm)+m), insize(dm));
                fns = 1;
            elseif(insize(dm)>1 && fns==1)
                idcell{dm} = ':';
            else
                idcell{dm} = 1;
            end
        end
    else

        for dm = 1:dms
            if(insize(dm)>1)
                idcell{dm} = mirroridx((-m(dm)+1):(insize(dm)+m(dm)), insize(dm));
            else
                idcell{dm} = 1;
            end
        end
    
    end
    
    sb.type = '()';
    sb.subs = idcell;
    
    y = subsref(x, sb);
    
end