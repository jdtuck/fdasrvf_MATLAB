function [q2best,Rbest,gamIbest,minE] = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag,rotFlag,closed,scale,lam,method)
% FIND_ROTATION_AND_SEED_UNIQUE find roation and seed of two curves
% -------------------------------------------------------------------------
% find roation and seed of closed curve
%
% Usage: [q2best,Rbest] = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag)
%
%
% Input:
% q1: matrix (n,T) defining T points on n dimensional SRVF
% q2: matrix (n,T) defining T points on n dimensional SRVF
% reparamFlag: flag to calculate reparametrization (default F)
% closed: flag if closed curve (default F)
% scale: scale curves to unit length
% lam: penalty (default 0.0)
% method: controls which optimization method (default="DP") options are
% Dynamic Programming ("DP") and Riemannian BFGS
% ("RBFGSM")
%
% Output:
% q2best: best aligned and rotated and seeded SRVF
% Rbest: best rotation matrix
% gamI: best reparameterization
% minE: minimum Energy

arguments
    q1
    q2
    reparamFlag = true;
    rotFlag = true;
    closed = false;
    scale = true;
    lam = 0.0;
    method = 'DP';
end

[n,T] = size(q1);

scl = 4;
minE = 1000;
if closed
    end_idx = floor((T)/scl);
else
    end_idx = 0;
end

for ctr = 0:end_idx
    if closed
        q2n = ShiftF(q2,scl*ctr);
    else
        q2n = q2;
    end
    
    if (rotFlag)
        [q2n,R] = Find_Best_Rotation(q1,q2n);
    else
        R = eye(n);
    end
    
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = optimum_reparam_curve(q2,q1,lam,method);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            q2new = warp_srvf_gamma(q2n,gamI,scale);
            if closed
                q2new = ProjectC(q2new);
            end
        else
            q2new = q2n;
            gamI = linspace(0,1,T);
        end
        
    else
        q2new  = q2n;
    end

    if scale
        tmp = InnerProd_Q(q1,q2new);
        if tmp > 1
            tmp = 1;
        elseif tmp < -1
            tmp = -1;
        end
        Ec = acos(tmp);
    else
        v = q1 - q2new;
        Ec = sqrt(InnerProd_Q(v,v));
    end
    
    if Ec < minE
        Rbest = R;
        q2best = q2new;
        gamIbest = gamI;
        minE = Ec;
    end
end
