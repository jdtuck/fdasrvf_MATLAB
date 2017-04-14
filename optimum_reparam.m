function gam = optimum_reparam(q1,q2,t,lambda,method,w,f1o,f2o)
addpath(genpath('DP'))
addpath(genpath('gropt'))
if nargin < 4
    lambda = 0.0;
    method = 'DP';
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
        [opt,swap,tmp,tmp1] = ElasticCurvesReparam(c1, c2, w, onlyDP,  ...
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
