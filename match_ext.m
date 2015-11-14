function [D,P, mpath] = match_ext(t1,ext1,d1,t2,ext2,d2)
    te1 = t1(ext1);
    te2 = t2(ext2);
    
    % We'll pad each sequence to start on a 'peak' and end on a 'valley'
    pad1 = [0,0];
    pad2 = [0,0];
    if d1==-1
        te1(2:end+1) = te1;
        te1(1) = te1(2);
        pad1(1) = 1;
    end
    if mod(length(te1),2)==1
        te1(end+1) = te1(end);
        pad1(2) = 1;
    end
    
    if d2==-1
        te2(2:end+1) = te2;
        te2(1) = te2(2);
        pad2(1) = 1;
    end
    if mod(length(te2),2)==1
        te2(end+1) = te2(end);
        pad2(2) = 1;
    end
    
    n1 = length(te1);
    n2 = length(te2);
    
    % initialize weight and path matrices
    D = zeros(n1,n2);
    P = zeros(n1,n2,2);
    
    for i=1:n1
        for j=1:n2
            if mod(i+j,2)==0
                for ib = (i-1):-2:(1+mod(i,2))
                    for jb = (j-1):-2:(1+mod(j,2))
                        icurr = ib+1:2:i;
                        jcurr = jb+1:2:j;
                        W = sqrt(sum(te1(icurr)-te1(icurr-1)))...
                            *sqrt(sum(te2(jcurr)-te2(jcurr-1)));
                        Dcurr = D(ib,jb) + W;
                        if Dcurr > D(i,j);
                            D(i,j) = Dcurr;
                            P(i,j,:) = [ib,jb];
                        end
                    end
                end
            end
        end
    end
    
    D = D(1+pad1(1):end-pad1(2), 1+pad2(1):end-pad2(2));
    P = P(1+pad1(1):end-pad1(2), 1+pad2(1):end-pad2(2), :);
    P(:,:,1) = P(:,:,1)-pad1(1);
    P(:,:,2) = P(:,:,2)-pad2(1);

    % retrieve best path
    
    if pad1(2)==pad2(2)
        mpath = size(D);
    elseif D(end-1,end) > D(end,end-1)
        mpath = size(D) - [1,0];
    else
        mpath = size(D) - [0,1];
    end
    prev_vert = reshape(P(mpath(1,1),mpath(1,2),:),1,2);
    while all(prev_vert>0);
        mpath = [prev_vert; mpath];
        prev_vert = reshape(P(mpath(1,1),mpath(1,2),:),1,2);
    end
    
end

