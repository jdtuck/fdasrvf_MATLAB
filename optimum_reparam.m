function gam = optimum_reparam(q1,q2,t,lambda,method,w,f1o,f2o)
% OPTIMUM_REPARAM Calculates Warping for two SRVFs
% -------------------------------------------------------------------------%
% This function aligns two SRSF functions using Dynamic Programming
%
% Usage:  gam = optimum_reparam(q1,q2,t)
%         gam = optimum_reparam(q1,q2,t,lambda)
%         gam = optimum_reparam(q1,q2,t,lambda,method)
%         gam = optimum_reparam(q1,q2,t,lambda,method,w,f1o,f2o)
%
% Input:
% q1: srsf of function 1
% q2: srsf of function 2
% t: sample points of function 2
% lambda: controls amount of warping (default = 0)
% method: controls which optimization method (default="DP") options are
% Dynamic Programming ("DP"), Coordinate Descent ("DP2"), and Riemannian BFGS
% ("RBFGS")
% w: controls LRBFGS (default = 0.01)
% f1o: initial value of f1, vector or scalar depending on q1, defaults to zero
% f2o: initial value of f2, vector or scalar depending on q1, defaults to zero
%
% Output:
% gam: warping function
if nargin < 4
    lambda = 0.0;
    method = 'DP1';
    w = 0.0;
    f1o = 0.0;
    f2o = 0.0;
elseif nargin < 5
    method = 'DP';
    w = 0.0;
    f1o = 0.0;
    f2o = 0.0;
elseif nargin < 6
    w = 0.0;
    f1o = 0.0;
    f2o = 0.0;
elseif nargin < 7
    f1o = 0.0;
    f2o = 0.0;
elseif nargin < 8
    f2o = 0.0;
end

q1 = q1/norm(q1);
q2 = q2/norm(q2);
c1 = srvf_to_f(q1,t,f1o);
c2 = srvf_to_f(q2,t,f2o);
rotated = 0;
isclosed = 0;
skipm = 0;
auto = 0;
switch upper(method)
    case 'DP'
        [G,T] = DynamicProgrammingQ2(q1',t',q2',t',t',t',lambda);
        gam0 = interp1(T,G,t);
    case 'DP1'
        [gam0] = DynamicProgrammingQ(q2',q1',lambda,0);
    case 'SIMUL'
        [s1,s2, g1,g2,~,~,~]  = simul_align(c1,c2);
        u = linspace(0,1,length(g1));
        tmin = min(t);
        tmax = max(t);
        t2 = t;
        t2 = (t2-tmin)/(tmax-tmin);
        gam0 = simul_gam(u,g1,g2,t2,s1,s2,t2);
    case 'DP2'
        onlyDP = 1;
        [opt,swap,~,~] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
            rotated, isclosed, skipm, 'LRBFGS', auto);
        gam0 = opt(1:end-2);
        
        if swap
            gam0 = invertGamma(gam0);
        end
    case 'RBFGS'
        onlyDP = 0;
        [opt,swap,fopts,~] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
            rotated, isclosed, skipm, 'LRBFGS', auto);
        
        
        if (fopts(1) == 1000)
            onlyDP = 1;
            [opt,swap,~,~] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
                rotated, isclosed, skipm, 'LRBFGS', auto);
        end
        
        gam0 = opt(1:end-2);
        
        if swap
            gam0 = invertGamma(gam0);
        end
    case 'RBFGSM'
        t1 = linspace(0,1,length(t));
        options.verbosity=0;
        options.maxiter=30;
        options.inittype='id';
        options.Tdp=t1;
        M = infspherefactory(t1);
        [q2,~]=initialize(q1,q2,M,options);
        
        
        [~,gam0,~, ~, ~]=rlbfgs(q1.',q2.',M,options);
        
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
        %         hdp=smooth(hdp,0.1,'loess'); % smooth 
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
            [q20,gam0,cost]=sgd(q1,q2,M,options);
    %         figure; plot(cost)
    %         keyboard;
            
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
