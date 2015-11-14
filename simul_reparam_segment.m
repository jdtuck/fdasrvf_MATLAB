function [gg1,gg2] = simul_reparam_segment(src, tgt, te1, te2)
    i1 = src(1)+1:2:tgt(1);
    i2 = src(2)+1:2:tgt(2);
    
    v1 = sum(te1(i1)-te1(i1-1));
    v2 = sum(te2(i2)-te2(i2-1));
    R = v2/v1;
    
    a1=src(1); a2=src(2);
    t1=te1(a1); t2=te2(a2);
    u1=0; u2=0;
    
    gg1=[]; gg2=[];
    
    while (a1<tgt(1)) && (a2<tgt(2))
        if a1==tgt(1)-1 && a2==tgt(2)-1
            a1=tgt(1); a2=tgt(2);
            gg1 = [gg1; te1(a1)];
            gg2 = [gg2; te2(a2)];
        else
            p1 = (u1 + te1(a1+1) - t1)/v1;
            p2 = (u2 + te2(a2+1) - t2)/v2;
            if p1<p2
                lam = t2 + R*(te1(a1+1)-t1);
                gg1 = [gg1; te1(a1+1); te1(a1+2)];
                gg2 = [gg2; lam; lam];
                u1 = u1 + te1(a1+1) - t1;
                u2 = u2 + lam - t2;
                t1 = te1(a1+2);
                t2 = lam;
                a1 = a1+2;
            else
                lam = t1 + (1/R)*(te2(a2+1)-t2);
                gg1 = [gg1; lam; lam];
                gg2 = [gg2; te2(a2+1); te2(a2+2)];
                u1 = u1 + lam - t1;
                u2 = u2 + te2(a2+1) - t2;
                t1 = lam;
                t2 = te2(a2+2);
                a2 = a2+2;
            end
        end
    end
end