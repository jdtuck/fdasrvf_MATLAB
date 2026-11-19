function [q2Opt,gammaOpt,cost,info,options] = rlbfgs(q1,q2,M,options)
% Riemannian limited memory BFGS solver for elastic registration, the MATLAB
% counterpart of c_rlbfgs (armadillo_cpp/rbfgs.hpp), which also handles curves
% in R^n.  It follows Huang, Gallivan, Srivastava and Absil, "Riemannian
% optimization for registration of curves in elastic shape analysis", for the
% warping alone (no rotation or shift).  The iterate is l = sqrt(gammadot) on
% the unit sphere of L^2([0,1],R), gamma = cumtrapz(l^2), and the cost is
%
%   ||q1 - q2(gamma) l||^2 + lam*penalty(l) + omega*B(l),
%
% where B is the paper's barrier, whose weight omega decays each iteration.
% Steps are taken along the exponential map and the stored BFGS pairs are moved
% to the new tangent space by parallel translation.
%
% function [q2Opt,gammaOpt,cost,info,options] = rlbfgs(q1,q2,M)
% function [q2Opt,gammaOpt,cost,info,options] = rlbfgs(q1,q2,M,options)
%
% The inputs q1 and q2 are the n-by-T square root velocity functions to be
% aligned; q2 is aligned to q1. M is the manifold factory for the unit
% hypersphere in L^2([0,1],R), on a uniform grid.
%
% The outputs 'gammaOpt' and 'cost' are the optimal warping and its cost
% without the barrier. q2Opt is q2 optimally warped to q1.
%
% The output 'info' is a struct-array with an element per iteration:
%   iter (integer)
%       The iteration number. The initial guess is 0.
%   cost (double)
%       The cost value, without the barrier.
%   gradnorm (double)
%       The (Riemannian) norm of the gradient.
%   time (double)
%       The total elapsed time in seconds to reach the corresponding cost.
%   stepsize (double)
%       The size of the step from the previous to the new iterate.
%   accepted (Boolean)
%       true if the step is accepted in the cautious update.
%   stop (integer)
%       Why the solver stopped, in the last element and 0 before: 1 no
%       acceptable step, 2 gradient norm, 3 iterations, 4 relative change of
%       the cost, 5 time.
%
% The options structure overwrites these defaults:
%
%   gam0 ([])
%       The warping to start from, sampled on M.t; empty starts from the
%       identity.
%   lam (0), penalty (0)
%       The penalty weight and type: 0 roughness, 1 l2gam, 2 l2psi, 3 geodesic.
%   maxiter (100)
%       The algorithm terminates if maxiter iterations were executed.
%   miniter (10), tolrelcost (1e-3), relcostcount (5)
%       The algorithm terminates, after at least miniter iterations, once the
%       relative change of the cost has stayed below tolrelcost for more than
%       relcostcount iterations.
%   tolgradnorm (1e-10)
%       The algorithm terminates if the norm of the gradient drops below this.
%   omega0 (1/100), omegarate (0.8)
%       The initial barrier weight and the factor it is cut by each iteration.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%   memory (30)
%       The number of previous iterations the program remembers.
%   strict_inc_func (@(t) 1e-4*t)
%       The Cautious step needs a real function that has value 0 at t = 0,
%       and is strictly increasing. See details in Wen Huang's paper
%       "A Riemannian BFGS Method without Differentiated Retraction for
%       Nonconvex Optimization Problems"
%   ls_max_steps (25)
%       The most cost evaluations of a line search.
%   plotevol (false)
%       This flag can be set true to plot the warped q2 at each iteration
%       along with a fixed q1.
%   verbosity (3)
%       Integer number used to tune the amount of output the algorithm
%       generates. The higher, the more output. 0 means silent.
%
% Please cite the Manopt paper as well as the research paper:
% @InBook{Huang2016,
%   title     = {A {R}iemannian {BFGS} Method for Nonconvex Optimization Problems},
%   author    = {Huang, W. and Absil, P.-A. and Gallivan, K.A.},
%   year      = {2016},
%   publisher = {Springer International Publishing},
%   editor    = {Karas{\"o}zen, B{\"u}lent and Manguo{\u{g}}lu, Murat and Tezer-Sezgin, M{\"u}nevver and G{\"o}ktepe, Serdar and U{\u{g}}ur, {\"O}m{\"u}r},
%   address   = {Cham},
%   booktitle = {Numerical Mathematics and Advanced Applications ENUMATH 2015},
%   pages     = {627--634},
%   doi       = {10.1007/978-3-319-39929-4_60}
% }

% This file is a modified version of rlbfgs.m in Manopt: www.manopt.org.
% Original author: Changshuo Liu, July 19, 2017.
% Contributors: Nicolas Boumal
% Modified by: Darshan Bryner, August 1, 2019.

% Local defaults for the program
localdefaults.gam0 = [];
localdefaults.lam = 0;
localdefaults.penalty = 0;
localdefaults.maxiter = 100;
localdefaults.miniter = 10;
localdefaults.tolrelcost = 1e-3;
localdefaults.relcostcount = 5;
localdefaults.tolgradnorm = 1e-10;
localdefaults.omega0 = 1/100;
localdefaults.omegarate = 0.8;
localdefaults.maxtime = inf;
localdefaults.memory = 30;
localdefaults.strict_inc_func = @(t) 1e-4*t;
localdefaults.ls_max_steps = 25;
localdefaults.plotevol = false;
localdefaults.verbosity = 3;

% Merge local defaults w/ user specified options, if any.
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);
options.memory = max(options.memory, 0);

timetic = tic();

% The problem: quantities that stay fixed through the iterations.  q1 and q2
% have unit L2 norm on M.t, so the cost and its gradient do not scale with the
% grid size, and lam is rescaled so that it weighs the penalty as it does for
% unit Euclidean norm SRVFs.
t = M.t;
T = M.T;
P.t = t;
P.T = T;
P.hb = 1/(T-1);
P.w = [P.hb/2, P.hb*ones(1,T-2), P.hb/2];
P.q1 = q1/sqrt(trapz(t,sum(q1.^2,1)));
P.q2 = q2/sqrt(trapz(t,sum(q2.^2,1)));
P.q2M = spline_second_derivs(P.q2,P.hb);
P.lam = options.lam*(norm(q1,'fro')/sqrt(trapz(t,sum(q1.^2,1))))* ...
    (norm(q2,'fro')/sqrt(trapz(t,sum(q2.^2,1))));
P.penalty = options.penalty;

l = ones(1,T);
if ~isempty(options.gam0)
    % sqrt(gammadot) from central differences; a small floor keeps l inside
    % the barrier
    gam = options.gam0(:).';
    gam = (gam-gam(1))/(gam(T)-gam(1));
    l = sqrt(max(fd_gradient(gam,P.hb),1e-6));
end
l = l/M.norm(l);

sHistory = cell(1,options.memory);
yHistory = cell(1,options.memory);
rhoHistory = zeros(1,options.memory);

j = 0;                    % number of stored pairs
k = 0;                    % total number of BFGS iterations
scaleFactor = 1;          % scaling of initial matrix, Barzilai-Borwein
stop = 0;
smallchange = 0;
failedsteps = 0;
omega = options.omega0;
stepsize = NaN;
accepted = NaN;

cur = costgrad(l,P);
f = cur.fixedcost + omega*cur.barrier;
g = riemannian_grad(l,cur,omega,P,M);
gnorm = M.norm(g);

stats = savestats();
info = stats;

if options.verbosity >= 2
    fprintf(' iter                   cost val            grad. norm\n');
end

while true

    if options.verbosity >= 2
        fprintf('%5d    %+.16e        %.8e\n',k,cur.fixedcost,gnorm);
    end

    if gnorm <= options.tolgradnorm
        stop = 2;
        break;
    end
    if k >= options.maxiter
        stop = 3;
        break;
    end
    if toc(timetic) >= options.maxtime
        stop = 5;
        break;
    end

    % compute BFGS direction; fall back to steepest descent if the
    % quasi-Newton direction is not a descent direction
    p = getDirection(M,g,sHistory,yHistory,rhoHistory,scaleFactor,j);
    df0 = M.inner(g,p);
    if df0 >= 0
        j = 0;
        scaleFactor = 1;
        p = -g;
        df0 = -gnorm^2;
    end

    % Armijo backtracking along the exponential map
    alpha = 1;
    lnew = M.exp(l,p,alpha);
    nxt = costgrad(lnew,P);
    fnew = nxt.fixedcost + omega*nxt.barrier;
    evals = 1;
    while any(lnew <= 0) || ~isfinite(fnew) || fnew > f + 1e-4*alpha*df0
        alpha = alpha*0.5;
        lnew = M.exp(l,p,alpha);
        nxt = costgrad(lnew,P);
        fnew = nxt.fixedcost + omega*nxt.barrier;
        evals = evals + 1;
        if evals >= options.ls_max_steps
            break;
        end
    end
    if any(lnew <= 0) || ~isfinite(fnew) || fnew > f
        % no acceptable step: restart the memory once, then stop
        failedsteps = failedsteps + 1;
        if failedsteps >= 2 || j == 0
            stop = 1;
            break;
        end
        j = 0;
        scaleFactor = 1;
        continue;
    end
    failedsteps = 0;
    stepsize = alpha*M.norm(p);

    % the gradient at the new point with this iteration's omega, so that the
    % BFGS pair is that of one cost function
    gnew = riemannian_grad(lnew,nxt,omega,P,M);

    % parallel translation of tangent vectors from l to lnew
    lsum = l + lnew;
    nsum = M.inner(lsum,lsum);
    transp = @(v) v - (2*M.inner(v,lnew)/nsum)*lsum;
    sk = transp(alpha*p);
    yk = gnew - transp(g);
    for i = 1:j
        sHistory{i} = transp(sHistory{i});
        yHistory{i} = transp(yHistory{i});
    end

    % cautious update
    inner_sk_yk = M.inner(sk,yk);
    inner_sk_sk = M.inner(sk,sk);
    accepted = inner_sk_sk > 0 && ...
        inner_sk_yk/inner_sk_sk >= options.strict_inc_func(gnorm);
    if accepted && options.memory > 0
        if j >= options.memory
            sHistory = sHistory([2:end,1]);
            yHistory = yHistory([2:end,1]);
            rhoHistory = rhoHistory([2:end,1]);
            j = options.memory - 1;
        end
        sHistory{j+1} = sk;
        yHistory{j+1} = yk;
        rhoHistory(j+1) = 1/inner_sk_yk;
        scaleFactor = inner_sk_yk/M.inner(yk,yk);
        j = j + 1;
    end

    % relative change of the cost, both values at this omega
    relchange = abs(f-fnew)/max(abs(f),1e-300);
    if relchange < options.tolrelcost
        smallchange = smallchange + 1;
    else
        smallchange = 0;
    end

    l = lnew;
    cur = nxt;
    k = k + 1;

    % the barrier fades as the iterations proceed
    omega = omega*options.omegarate;
    f = cur.fixedcost + omega*cur.barrier;
    g = riemannian_grad(l,cur,omega,P,M);
    gnorm = M.norm(g);

    if options.plotevol
        figure(100); clf; plot(t,P.q1,'b',t,warp_q(l,P),'r');
        pause;
    end

    stats = savestats();
    info(k+1) = stats;

    if k >= options.miniter && smallchange > options.relcostcount
        stop = 4;
        break;
    end
end

% Housekeeping before we return
info(end).stop = stop;
gam = cumtrapz(t,l.^2);
gammaOpt = (gam-min(gam))/(max(gam)-min(gam));
q2Opt = warp_q(l,P);
cost = cur.fixedcost;

if options.verbosity >= 1
    fprintf('Stopped with reason %d after %d iterations in %f [s]\n', ...
        stop,k,toc(timetic));
end

% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = k;
        stats.cost = cur.fixedcost;
        stats.gradnorm = gnorm;
        stats.time = toc(timetic);
        stats.stepsize = stepsize;
        stats.accepted = accepted;
        stats.stop = 0;
    end

end


function dir = getDirection(M,g,sHistory,yHistory,rhoHistory,scaleFactor,j)
% BFGS step, see Wen's paper for details. This function takes in a tangent
% vector g, and applies an approximate inverse Hessian to it to get Pg.  Then,
% -Pg is returned.

q = g;
inner_s_q = zeros(1,j);

for i = j:-1:1
    inner_s_q(i) = rhoHistory(i)*M.inner(sHistory{i},q);
    q = q - inner_s_q(i)*yHistory{i};
end

r = scaleFactor*q;

for i = 1:j
    om = rhoHistory(i)*M.inner(yHistory{i},r);
    r = r + (inner_s_q(i)-om)*sHistory{i};
end

dir = -r;

end

function c = costgrad(l,P)
% Cost parts and their partial derivatives with respect to each sample of l,
% differentiating the discretised cost exactly: gamma is cumtrapz(l^2), q2 is
% evaluated on its cubic spline and the integrals are trapezoidal sums.  The
% part that stays fixed (data and penalty) is split from the barrier, whose
% weight changes each iteration.

w = P.w;
gam = cumtrapz(P.t,l.^2);
[Q,Qdot] = spline_eval(P.q2,P.q2M,gam,P.hb);
r = P.q1 - Q.*l;
c.fixedcost = sum(w.*sum(r.^2,1));
% r depends on l directly and through gamma
a = -2*(w.*sum(r.*Qdot,1).*l);
c.fixedgrad = -2*(w.*sum(r.*Q,1)) + 2*(l.*cumtrapz_adjoint(a,P.hb));

if P.lam ~= 0
    c.fixedcost = c.fixedcost + P.lam*penalty_cost(l,P);
    c.fixedgrad = c.fixedgrad + P.lam*penalty_grad(l,P);
end

% B = trapz((l^2 + 1/l^2) sqrt(1 + l^4)), eq. (15) of the paper; it grows where
% gamma is nearly flat or vertical and is symmetric under inverting gamma
l2 = l.^2;
rt = sqrt(1 + l.^4);
u = l2 + 1./l2;
c.barrier = sum(w.*(u.*rt));
c.barriergrad = w.*((2*l - 2./l.^3).*rt + u.*(2*l.^3./rt));

end

function g = riemannian_grad(l,c,omega,P,M)
% Riemannian gradient under the L2 metric M.inner uses: the partial derivatives
% divided by the quadrature weights, projected onto the tangent space at l

v = (c.fixedgrad + omega*c.barriergrad)./P.w;
g = v - M.inner(v,l)*l;

end

function qw = warp_q(l,P)
% q2 warped by the warping whose SRVF is l

Q = spline_eval(P.q2,P.q2M,cumtrapz(P.t,l.^2),P.hb);
qw = Q.*l;

end

function pen = penalty_cost(h,P)
% The penalty at the warping whose SRVF is h: 0 roughness, 1 l2gam, 2 l2psi,
% 3 geodesic

w = P.w;
switch P.penalty
    case 0
        pen = sum(w.*fd_gradient(h.^2,P.hb).^2);
    case 1
        pen = sum(w.*(h.^2-1).^2);
    case 2
        pen = sum(w.*(h-1).^2);
    case 3
        pen = acos(min(max(sum(w.*h),-1),1))^2;
end

end

function dP = penalty_grad(h,P)
% Partial derivatives of penalty_cost with respect to each sample of h

w = P.w;
switch P.penalty
    case 0
        hdot = fd_gradient(h.^2,P.hb);
        dP = 2*(h.*gradient_adjoint(2*(w.*hdot),P.hb));
    case 1
        dP = 4*(h.*(h.^2-1)).*w;
    case 2
        dP = 2*(h-1).*w;
    case 3
        c = min(max(sum(w.*h),-1),1);
        % acos(c)/sqrt(1-c^2) tends to 1 as c tends to 1, where both factors
        % vanish
        scale = 1;
        if 1 - c^2 > 1e-12
            scale = acos(c)/sqrt(1-c^2);
        end
        dP = -2*scale*w;
end

end

function g = fd_gradient(f,hb)
% central differences inside, one sided at the ends

g = [f(2)-f(1), (f(3:end)-f(1:end-2))/2, f(end)-f(end-1)]/hb;

end

function out = gradient_adjoint(y,hb)
% the transpose of the finite difference operator fd_gradient applies

T = numel(y);
out = zeros(1,T);
out(2) = out(2) + y(1)/hb;
out(1) = out(1) - y(1)/hb;
out(T) = out(T) + y(T)/hb;
out(T-1) = out(T-1) - y(T)/hb;
out(3:T) = out(3:T) + y(2:T-1)/(2*hb);
out(1:T-2) = out(1:T-2) - y(2:T-1)/(2*hb);

end

function out = cumtrapz_adjoint(b,hb)
% the transpose of the running trapezoidal sum on a uniform grid

T = numel(b);
acc = cumsum(b(T:-1:2));
acc = acc(end:-1:1)*hb/2;
out = zeros(1,T);
out(2:T) = acc;
out(1:T-1) = out(1:T-1) + acc;

end

function Mq = spline_second_derivs(y,hb)
% second derivatives of the natural cubic spline interpolating the rows of y

[n,T] = size(y);
Mq = zeros(n,T);
if T < 3
    return
end
e = ones(T-2,1);
A = spdiags([e 4*e e],-1:1,T-2,T-2);
rhs = 6*(y(:,3:T)-2*y(:,2:T-1)+y(:,1:T-2))/hb^2;
Mq(:,2:T-1) = (A\rhs.').';

end

function [val,dval] = spline_eval(y,Mq,x,hb)
% value and derivative at x of the natural cubic spline of y, whose second
% derivatives are Mq

T = size(y,2);
xi = min(max(x,0),1);
k = min(floor(xi/hb),T-2);
A = (k+1)*hb - xi;   % distance to the right node
B = xi - k*hb;       % distance to the left node
yk = y(:,k+1);
yk1 = y(:,k+2);
Mk = Mq(:,k+1);
Mk1 = Mq(:,k+2);
cA = yk/hb - Mk*hb/6;
cB = yk1/hb - Mk1*hb/6;
val = Mk.*A.^3/(6*hb) + Mk1.*B.^3/(6*hb) + cA.*A + cB.*B;
dval = -Mk.*A.^2/(2*hb) + Mk1.*B.^2/(2*hb) - cA + cB;

end

function opts = mergeOptions(opts_sub, opts_master)
% Merges two options structures with one having precedence over the other:
% whenever a field is present in both, it is the value in opts_master that is
% kept.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.

if isempty(opts_sub)
    opts_sub = struct();
end
if isempty(opts_master)
    opts_master = struct();
end

opts = opts_sub;
fields = fieldnames(opts_master);
for i = 1 : length(fields)
    opts.(fields{i}) = opts_master.(fields{i});
end

end
