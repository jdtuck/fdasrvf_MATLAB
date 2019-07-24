function [c,ctilde] = basis_tangent_id(n,t)

T=length(t);
N=2*n;
c=zeros(N,T);
ctilde=zeros(N,T);
j=0;
for i=1:n
    j=j+1;
    c(j,:)=2*sin(2*pi*i*t)/sqrt(2);
    ctilde(j,:)=cumtrapz(t,c(j,:));
    j=j+1;
    c(j,:)=2*cos(2*pi*i*t)/sqrt(2);
    ctilde(j,:)=cumtrapz(t,c(j,:));
end