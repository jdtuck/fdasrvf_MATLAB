function [q2opt,gammaOpt,cost] = sgd(q1,q2,M,options)

localdefaults.maxepochs=10;
localdefaults.plotevol=false;
localdefaults.stepsize=0.001;
localdefaults.maxfreq=10;

% Merge local defaults w/ user specified options, if any.
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults,options);

t=M.t;
T=M.T;
[p,~]=size(q1);
hid=ones(1,T);
n=options.maxfreq;
% n=round(T/2)-1;
[c,ctilde]=basis_tangent_id(n,t);
[N,~]=size(c);
alpha0=options.stepsize;
K=options.maxepochs;
htilde=zeros(N*K+1,T);
htilde(1,:)=hid;
q2tilde=zeros(N*K+1,T);
q2tilde(1,:)=q2;
cost=zeros(K*N+1,1);
cost(1)=alignment_cost(hid,q1,q2tilde(1,:),M);
iter=1;
if options.plotevol
    figure(100); clf; plot(t,q1,t,q2tilde(iter,:));
%         figure(101); clf; plot(M.t,htilde(iter,:));
%             pause;
end
% keyboard;
for k=1:K
    idx=randperm(N);
    for j=1:N
        i=idx(j);
        q2tildedot=gradient(q2tilde(iter,:),1/(T-1));
        vi=c(i,:)*innerProdL2(q1-q2tilde(iter,:),2*q2tildedot.*repmat(ctilde(i,:),p,1)+q2tilde(iter,:).*repmat(c(i,:),p,1),t);
        alpha=alpha0;
        hi=M.exp(hid,vi,alpha);
        while sum(hi<0)>0
            alpha=alpha/2;
            hi=M.exp(hid,vi,alpha);
        end
        htilde(iter+1,:)=group_action_SRVF(htilde(iter,:),hi,M);
        q2tilde(iter+1,:)=group_action_SRVF(q2tilde(iter,:),hi,M);
        iter=iter+1;
        cost(iter)=alignment_cost(hid,q1,q2tilde(iter,:),M);
        if options.plotevol
            figure(100); clf; plot(t,q1,t,q2tilde(iter,:));
    %         figure(101); clf; plot(M.t,htilde(iter,:));
%             pause;
        end
    end
end
[~,idx]=min(cost);
q2opt=q2tilde(idx,:);
gammaOpt=cumtrapz(t,htilde(idx,:).^2);
gammaOpt=gammaOpt/gammaOpt(end);
        
