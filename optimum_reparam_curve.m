function gam = optimum_reparam_curve(q1,q2,lam,method)
% OPTIMUM_REPARAM_CURVE Calculates Warping for two SRVFs
% -------------------------------------------------------------------------%
% This function aligns two SRVF functions using Dynamic Programming
%
% Usage:  gam = optimum_reparam_curve(q1,q2)
%         gam = optimum_reparam_curve(q1,q2,method)
%
% Input:
% q1: srvf of function 1
% q2: srvf of function 2
% lam: penalty (default 0.0)
% method: controls which optimization method (default="DP") options are
% Dynamic Programming ("DP") and Riemannian BFGS
% ("RBFGSM")
%
% Output:
% gam: warping function
if nargin < 3
    lam = 0.0;
    method = 'DP';
elseif nargin < 4
    method = 'DP';
end

q1 = q1/sqrt(InnerProd_Q(q1,q1));
q2 = q2/sqrt(InnerProd_Q(q2,q2));

switch upper(method)
    case 'DP'
        gam0 = DynamicProgrammingQ(q2,q1,lam,0);
    case 'RBFGSM'
        t1 = linspace(0,1,size(q1,2));
        options.verbosity=0;
        options.maxiter=30;
        options.inittype='id';
        options.Tdp=t1;
        M = infspherefactory(t1);
        [q2,~]=initialize(q1,q2,M,options);
        
        [~,gam0,~, ~, ~]=rlbfgs(q1,q2,M,options);
        
    otherwise
        error('Invalid Method')
end

gam = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end


% helper functions
function [q20,gam0] = initialize(q1,q2,M,options)

    t=M.t;
    T=M.T;
    
    switch options.inittype
        
        % Initialize as identity
        case 'id' 
            gam0=t;
            q20=q2;
            
        % Initialize by randomly warping q2
        case 'rand' 
            hid=ones(1,T);
            [b,~]=basis_tangent_id(1,t);
            [N,~]=size(b);
            Ninit=options.Ninit;
            htry=zeros(Ninit,T);
            costtry=zeros(Ninit,1);
            for i=1:Ninit
                coeffs=0.2*randn(1,N);
                v=coeffs*b;
                htry(i,:)=M.exp(hid,v,1);
                while sum(htry(i,:)<=0)>0
                    coeffs=0.2*randn(1,N);
                    v=coeffs*b;
                    htry(i,:)=M.exp(hid,v,1);
                end
                costtry(i)=alignment_cost(htry(i,:),q1,q2,M);
            end
            [~,idx]=min(costtry);
            h0=htry(idx,:);
            [q20,gam0]=group_action_SRVF(q2,h0,M);
            
        % Initialize from dynamic programming solution 
        case 'dp'
            
            Tdp=options.Tdp;
        
            % Resample
            if Tdp~=T
                tdp=linspace(0,1,Tdp);
                q1dn=spline(M.t,q1,tdp);
                q2dn=spline(M.t,q2,tdp);
                gamdp=DynamicProgrammingQ(q2dn,q1dn,0,1);
                hdp=sqrt(gradient(gamdp,1/(Tdp-1)));
                h0=spline(tdp,hdp,t); % resample 
                [q20,gam0]=group_action_SRVF(q2,h0,M);
            else
                gam0=DynamicProgrammingQ(q2,q1,0,1);
                [n,~]=size(q2);
                q20=spline(t,q2,gam0).*repmat(sqrt(gradient(gam0,1/(T-1))),n,1);
            end
           
        case 'sgd'
            
            options.plotevol=false;
            options.stepsize=0.001;
            [q20,gam0,~]=sgd(q1,q2,M,options);

            
    end
    
    
end    

function f = alignment_cost(h,q1,q2k,M)

    % Evaluate the cost function f = ||q1 - ((q2,hk),h)||^2.
    % h=sqrt{\dot{\gamma}} is a sequential update of cumulative warping hk
    
    t=M.t;
    q2new=group_action_SRVF(q2k,h,M);
    f=normL2(q1-q2new,t)^2; % Cost
end

function val = normL2(f,t)

    val=sqrt(innerProdL2(f,f,t));
end
