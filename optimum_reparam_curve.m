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
arguments
    q1
    q2
    lam = 0.0;
    method = 'DP';
end

q1 = q1/sqrt(InnerProd_Q(q1,q1));
q2 = q2/sqrt(InnerProd_Q(q2,q2));

switch upper(method)
    case 'DP'
        gam0 = DynamicProgrammingQ(q2,q1,lam,0);
    case 'RBFGSM'
        % As in Huang et al., "Riemannian optimization for registration of
        % curves in elastic shape analysis": start from a DP warp on a coarse
        % grid of at most 51 points and refine it on q1 and q2 resampled
        % linearly 3 times finer, where the trapezoid rule resolves q2 under
        % a steep warp
        M = size(q1,2);
        t1 = linspace(0,1,M);
        tc = linspace(0,1,min(M,51));
        [G,T] = DynamicProgrammingQ2(spline(t1,q1,tc),tc,spline(t1,q2,tc),tc,tc,tc,0,10,1);
        scl = 3;
        tf = linspace(0,1,(M-1)*scl+1);
        options.verbosity = 0;
        options.maxiter = 30;
        if ~isequal(G,T)
            % otherwise DP chose the identity, where the solver starts anyway
            options.gam0 = interp1(T,G,tf);
        end
        [~,gamf] = rlbfgs(interp1(t1,q1',tf')',interp1(t1,q2',tf')',infspherefactory(tf),options);
        % every scl-th point of tf is a point of t1
        gam0 = gamf(1:scl:end);
    otherwise
        error('Invalid Method')
end

gam = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end
