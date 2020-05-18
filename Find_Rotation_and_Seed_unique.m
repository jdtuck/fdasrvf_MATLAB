function [q2best,Rbest,gamIbest] = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag,closed)
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
%
% Output:
% q2best: best aligned and rotated and seeded SRVF
% Rbest: best rotation matrix
% gamI: best reparameterization

if nargin < 3
    reparamFlag = false;
    closed = false;
end

[~,T] = size(q1);

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
    [q2n,R] = Find_Best_Rotation(q1,q2n);
    
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1/sqrt(InnerProd_Q(q1,q1)),q2n/sqrt(InnerProd_Q(q2n,q2n)),0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            p2n = q_to_curve(q2n);
            p2new = warp_curve_gamma(p2n,gamI);
            q2new = curve_to_q(p2new);
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
    Ec = acos(InnerProd_Q(q1,q2new));
    if Ec < minE
        Rbest = R;
        q2best = q2new;
        gamIbest = gamI;
        minE = Ec;
    end
end
