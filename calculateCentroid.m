function centroid = calculateCentroid(beta)
% CALCULATECENTROID Calculate centroid of curve
% -------------------------------------------------------------------------
% Calculate centroid of curve
% 
% Usage: centroid = calculateCentroid(beta)
%%
% Input:
% beta: matrix (n,T) defining T points on n dimensional curve
% 
% Output:
% centroid: centroid of curve

[n,T]=size(beta);
betadot=gradient(beta,1/(T-1));
normbetadot=zeros(1,T);
integrand=zeros(n,T);
for i=1:T
    normbetadot(i)=norm(betadot(:,i));
    integrand(:,i)=beta(:,i)*normbetadot(i);
end
scale=trapz(linspace(0,1,T),normbetadot);
centroid=trapz(linspace(0,1,T),integrand,2)/scale;