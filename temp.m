T=100; % Number of samples
n=8; % Number of basis elements / 2
lambda=0; % penalty coeff
warping_amount=0.2; % smaller number generates random warpings closer to identity
Ntrain=50; % Number of warpings to generate per signal (class)

t=linspace(0,1,T);
qgam_id=ones(1,T);
M=infspherefactory(t);

%% Create original signals

forig(1,:)=exp(-3*t).*sin(3*pi*t);
fdot=gradient(forig(1,:),1/(T-1));
forig(1,:)=forig(1,:)/trapz(t,abs(fdot));
[Nclass,~]=size(forig);

% Create basis for generating random warpings
[b,~]=basis_tangent_id(2,t);
[Nb,~]=size(b);

[Nclass,~]=size(forig);

% Create basis for generating random warpings
[b,~]=basis_tangent_id(2,t);
[Nb,~]=size(b);

% Generate random warpings
K=Ntrain*Nclass;
c=zeros(K,Nb);
v=zeros(K,T);
f=zeros(K,T);
x=zeros(K,T);
xS=zeros(K,T);
k=0;
for i=1:Nclass
    
    for j=1:Ntrain
        k=k+1;
        c(k,:)=warping_amount*randn(1,Nb);
        v(k,:)=c(k,:)*b;
        qgamtemp=M.exp(qgam_id,v(k,:));
        gamtrue=cumtrapz(t,qgamtemp.^2);
        gamtrue=gamtrue/gamtrue(end);
        c(k,:)=warping_amount*randn(1,Nb);
        v(k,:)=c(k,:)*b;
        f(k,:)=spline(t,forig(i,:),gamtrue);
        f(k,:)=f(k,:)/trapz(t,abs(gradient(f(k,:),1/(T-1))));
        x(k,:)=15*rand*f(k,:);
        xS(k,:)=x(k,:)/norm(x(k,:));
    end
    
end

%% Alignment

obj = fdawarp(xS',t);
obj1 = obj;
tic
obj = obj.time_warping();
toc

option.parallel = 0;
option.closepool = 0;
option.smooth = 0;
option.sparam = 25;
option.method = 'RBFGS';
option.w = 0.0;
option.MaxItr = 20;

tic
obj1 = obj1.time_warping(0,option);
toc