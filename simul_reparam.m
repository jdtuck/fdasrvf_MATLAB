function [g1, g2] = simul_reparam(te1, te2, mpath)
    g1(1) = 0;
    g2(1) = 0;
    
    if mpath(1,2) == 2
        g1 = [g1; 0];
        g2 = [g2; te2(2)];
    elseif mpath(1,1) == 2
        g1 = [g1; te1(2)];
        g2 = [g2; 0];
    end
    
    m = length(mpath);
    for i=1:m-1
        [gg1,gg2] = simul_reparam_segment(mpath(i,:), mpath(i+1,:), te1, te2);
        
        g1 = [g1; gg1];
        g2 = [g2; gg2];
    end
    
    n1 = length(te1);
    n2 = length(te2);
    if (mpath(end,1) == n1-1) || (mpath(end,2) == n2-1)
        g1 = [g1; 1];
        g2 = [g2; 1];
    end
end