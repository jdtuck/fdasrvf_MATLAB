function [v,dist,gamI] = inverseExp_Coord(beta1,beta2)
% Calculate the inverse exponential to obtain a shooting vector from beta1 to
% [beta2] in shape space of open curves

T=length(beta1);

beta1=beta1-repmat(calculateCentroid(beta1),1,T);
beta2=beta2-repmat(calculateCentroid(beta2),1,T);

q1=curve_to_q(beta1);

% Iteratively optimize over SO(n) x Gamma
for i=1:1
    % Optimize over SO(n) 
    [beta2,~]=Find_Rotation_and_Seed_Coord(beta1,beta2);
    q2=curve_to_q(beta2);

    % Optimize over Gamma
    gam = DynamicProgrammingQ(q1,q2,0,1);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    beta2 = warp_curve_gamma(beta2,gamI);
end

% Optimize over SO(n)
[beta2,~]=Find_Rotation_and_Seed_Coord(beta1,beta2);
q2n=curve_to_q(beta2);

% Compute geodesic distance
q1dotq2=InnerProd_Q(q1,q2n);
dist=acos(q1dotq2); 

% Compute shooting vector
if q1dotq2>1
    q1dotq2=1;
end
u=q2n-q1dotq2*q1;
normu=sqrt(InnerProd_Q(u,u));
if normu>10^-4
    v=u*acos(q1dotq2)/normu;
else
    v=zeros(size(beta1,1),T);
end
