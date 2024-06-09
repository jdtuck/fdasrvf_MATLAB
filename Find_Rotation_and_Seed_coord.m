function [beta2best,q2best,Rbest,gamIbest] = Find_Rotation_and_Seed_coord(beta1,beta2,reparamFlag,rotation,closed,method)
% FIND_ROTATION_AND_SEED_COORD find roation and seed of two curves
% -------------------------------------------------------------------------
% find roation and seed of curves
%
% Usage: [q2best,Rbest] = Find_Rotation_and_Seed_coord(q1,q2,reparamFlag)
%
%
% Input:
% beta1: matrix (n,T) defining T points on n dimensional curve
% beta2: matrix (n,T) defining T points on n dimensional curve
% reparamFlag: flag to calculate reparametrization (default T)
% rotation: flag to compute optimal rotation (default T)
% closed: flag if closed curve (default F)
% method: optimzation method (default "DP")
%
% Output:
% q2best: best aligned and rotated and seeded SRVF
% Rbest: best rotation matrix
% gamI: best reparameterization

if nargin < 3
    reparamFlag = true;
    rotation = true;
    closed = false;
    method = 'DP';
elseif nargin < 4
    rotation = true;
    closed = false;
    method = 'DP';
elseif nargin < 5
    closed = false;
    method = 'DP';
elseif nargin < 6
    method = 'DP';
end

[n,T] = size(beta1);
q1 = curve_to_q(beta1);

scl = 4;
minE = 1000;
if closed
    end_idx = floor((T)/scl);
else
    end_idx = 0;
end

for ctr = 0:end_idx
    if closed
        beta2n = ShiftF(beta2,scl*ctr);
    else
        beta2n = beta2;
    end
    
    if (rotation)
        [beta2n,R] = Find_Best_Rotation(beta1,beta2n);
    else
        R = eye(n);
    end
    q2n = curve_to_q(beta2n);
    
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = optimum_reparam_curve(q2n,q1,0,method);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            beta2n = warp_curve_gamma(beta2n,gamI);
            q2n = curve_to_q(beta2n);
            if closed
                q2n = ProjectC(q2n);
            end
            q2new = q2n;
            beta2new = beta2n;
        else
            q2new = q2n;
            beta2new = beta2n;
            gamI = linspace(0,1,T);
        end
        
    else
        q2new  = q2n;
        beta2new = beta2n;
    end
    Ec = acos(InnerProd_Q(q1,q2new));
    if Ec < minE
        Rbest = R;
        q2best = q2new;
        beta2best = beta2new;
        gamIbest = gamI;
        minE = Ec;
    end
end
