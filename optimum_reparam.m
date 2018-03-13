function gam = optimum_reparam(q1,q2,t,lambda,method,w,f1o,f2o)
% Align two functions
%
% This function aligns two SRSF functions using Dynamic Programming
%
% @param q1 srsf of function 1
% @param q2 srsf of function 2
% @param t sample points of function 2
% @param lambda controls amount of warping (default = 0)
% @param method controls which optimization method (default="DP") options are
% Dynamic Programming ("DP"), Coordinate Descent ("DP2"), and Riemannian BFGS
% ("RBFGS")
% @param w controls LRBFGS (default = 0.01)
% @param f1o initial value of f1, vector or scalar depending on q1, defaults to zero
% @param f2o initial value of f2, vector or scalar depending on q1, defaults to zero
% @return gam warping function
addpath(genpath('DP'))
addpath(genpath('gropt'))
if nargin < 4
    lambda = 0.0;
    method = 'DP1';
    w = 0.0;
    f1o = 0.0;
    f2o = 0.0;
elseif nargin < 5
    method = 'DP1';
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
auto = 4;
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
            rotated, isclosed, skipm, 'DP', auto);
        gam0 = opt(1:end-2);

        if swap
            gam0 = invertGamma(gam0);
        end
    case 'RBFGS'
        onlyDP = 0;
        [opt,swap,fopts,~] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
            rotated, isclosed, skipm, method, auto);


        if (fopts(1) == 1000)
            onlyDP = 1;
            [opt,swap,~,~] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
                rotated, isclosed, skipm, method, auto);
        end

        gam0 = opt(1:end-2);

        if swap
            gam0 = invertGamma(gam0);
        end
    otherwise
        error('Invalid Method')
end

gam = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
