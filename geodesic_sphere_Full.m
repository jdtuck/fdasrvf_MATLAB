function X = geodesic_sphere_Full(q1,q2,stp,closed)
% GEODESIC_SPHERE_FULL Calculates geodesic on sphere
% -------------------------------------------------------------------------
% This function calculates a geodesic on the sphere at the
% starting point x_init in the direction g
% 
% Usage: X = geodesic_sphere_Full(q1,q2,stp,closed)
%
%
% Input:
% q1: matrix (n,T) defining T points on n dimensional SRVF1
% q2: matrix (n,T) defining T points on n dimensional SRVF2
% stp: step size
% closed: if closed curve
% 
% Output:
% X: matrix of geodsic samples

theta = acos(InnerProd_Q(q1,q2));
[n,T] = size(q1);
X = zeros(n,T,stp+1);
if theta > 0.0001
    for t=1:stp+1
        tau = (t-1)/stp;
        X(:,:,t) = (sin((1-tau)*theta)*q1 + sin(tau*theta)*q2)/sin(theta);
        if closed
            X(:,:,t) = ProjectC(X(:,:,t));
        end
        
    end
else
    for t=1:stp+1
        X(:,:,t) = q1;
    end
end