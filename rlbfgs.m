function [q2Opt,gammaOpt,cost,info,options] = rlbfgs(q1,q2,M,options)
% Riemannian limited memory BFGS solver for elastic function registration.
% The solver is designed to operate on the positive orthant of the unit
% hypersphere in L^2([0,1],R). The set of all functions
% h=\sqrt{\dot{\gamma}}, where \gamma is a diffeomorphism, is that
% manifold.
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below.
%
% function [q2new,gammaOpt,cost,info,options] = alignment_rlbfgs(q1,q2,c,ctilde,M,lambda)
% function [q2new,gammaOpt,cost,info,options] = alignment_rlbfgs(q1,q2,c,ctilde,M,lambda,options)
%
% The inputs q1 and q2 are the square root velocity functions of curves in
% R^n to be aligned. Here, q2 will be aligned to q1. M is the manifold
% factory for the unit hypersphere in L^2([0,1],R).
%
%
% The two outputs 'gammaOpt' and 'cost' are the optimal diffeomorphism
% and its cost, respectively. q2new is the optimally warped q2 to q1.
%
% The output 'info' is a struct-array which contains information about the
% iterations:
%   iter (integer)
%       The iteration number. The initial guess is 0.
%	cost (double)
%       The corresponding cost value.
%	gradnorm (double)
%       The (Riemannian) norm of the gradient.
%	time (double)
%       The total elapsed time in seconds to reach the corresponding cost.
%	stepsize (double)
%       The size of the step from the previous to the new iterate.
%   accepted (Boolean)
%       true if step is accepted in the cautious update. 0 otherwise.
%   And possibly additional information logged by options.statsfun.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached at each iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tolgradnorm (1e-3)
%       The algorithm terminates if the norm of the gradient drops below
%       this. For well-scaled problems, a rule of thumb is that you can
%       expect to reduce the gradient norm by 8 orders of magnitude
%       (sqrt(eps)) compared to the gradient norm at a "typical" point (a
%       rough initial iterate for example). Further decrease is sometimes
%       possible, but inexact floating point arithmetic will eventually
%       limit the final accuracy. If tolgradnorm is set too low, the
%       algorithm may end up iterating forever (or at least until another
%       stopping criterion triggers).
%   maxiter (100)
%       The algorithm terminates if maxiter iterations were executed.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%   minstepsize (1e-10)
%     The minimum norm of the tangent vector that points from the current
%     point to the next point. If the norm is less than minstepsize, the
%     program will terminate.
%   memory (30)
%     The number of previous iterations the program remembers. This is used
%     to approximate the inverse Hessian at the current point. Because of
%     difficulty of maintaining a representation of operators in terms of
%     coordinates, a recursive method is used. The number of steps in the
%     recursion is at most options.memory. This parameter can take any
%     integer value >= 0, or Inf, which is taken to be options.maxiter. If
%     options.maxiter has value Inf, then it will take value 10000 and a
%     warning will be displayed.
%   strict_inc_func (@(t) 1e-4*t)
%     The Cautious step needs a real function that has value 0 at t = 0,
%     and  is strictly increasing. See details in Wen Huang's paper
%     "A Riemannian BFGS Method without Differentiated Retraction for
%     Nonconvex Optimization Problems"
%   inittype ('id')
%     The type of initialization to be input into the function
%     initialize.m. So far, options are 'id' the identity element, 'rand' a
%     random initialization generated in some optimal sense, and 'dp' the
%     solution to a downsampled version of the optimization problem using
%     dynamic programming.
%   plotevol (false)
%     This flag can be set true to plot the cumulative warped q2 at each
%     iteration along with a fixed q1.
%   statsfun (none)
%       Function handle to a function that will be called after each
%       iteration to provide the opportunity to log additional statistics.
%       They will be returned in the info struct. See the generic Manopt
%       documentation about solvers for further information. statsfun is
%       called with the point x that was reached last.
%   stopfun (none)
%       Function handle to a function that will be called at each iteration
%       to provide the opportunity to specify additional stopping criteria.
%       See the generic Manopt documentation about solvers for further
%       information.
%   verbosity (2)
%       Integer number used to tune the amount of output the algorithm
%       generates during execution (mostly as text in the command window).
%       The higher, the more output. 0 means silent. 3 and above includes a
%       display of the options structure at the beginning of the execution.
%   debug (false)
%       Set to true to allow the algorithm to perform additional
%       computations for debugging purposes. If a debugging test fails, you
%       will be informed of it, usually via the command window. Be aware
%       that these additional computations appear in the algorithm timings
%       too, and may interfere with operations such as counting the number
%       of cost evaluations, etc. (the debug calls get storedb too).
%
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
%
% We point out that, at the moment, this implementation of RLBFGS can be
% slower than the implementation in ROPTLIB by Wen Huang et al. referenced
% above. For the purpose of comparing to their work, please use their
% implementation.
%


% This file is a modified version of rlbfgs.m in Manopt: www.manopt.org.
% Original author: Changshuo Liu, July 19, 2017.
% Contributors: Nicolas Boumal
% Modified by: Darshan Bryner, August 1, 2019.

% Local defaults for the program
localdefaults.minstepsize = 1e-10;
localdefaults.maxiter = 100;
localdefaults.tolgradnorm = 1e-3;
localdefaults.memory = 30;
localdefaults.strict_inc_func = @(t) 1e-4*t;
localdefaults.ls_max_steps = 25;
localdefaults.plotevol = false;
localdefaults.verbosity = 3;
localdefaults.maxtime = inf;

% Merge local defaults w/ user specified options, if any.
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

% To make sure memory in range [0, Inf)
options.memory = max(options.memory, 0);
if options.memory == Inf
    if isinf(options.maxiter)
        options.memory = 10000;
        warning('rlbfgs:memory', ['options.memory and options.maxiter' ...
            ' are both Inf; options.memory has been changed to 10000.']);
    else
        options.memory = options.maxiter;
    end
end

timetic = tic();

% __________Initialization of variables______________

% Initialize sequential update variables
htilde = ones(1,M.T);
q2tilde = q2;

% Number of iterations since the last restart
j = 0;

% Total number of BFGS iterations
k = 0;

% This cell stores step vectors which point from h_id to h_{k+1} for k
% indexing the last iterations, capped at options.memory.
% That is, it stores up to options.memory of the most recent step
% vectors in the tangent space of identity
sHistory = cell(1, options.memory);

% This cell stores the differences for latest k's of the gradient at
% time k+1 and the gradient at time k. The memory is also capped at
% options.memory.
yHistory = cell(1, options.memory);

% rhoHistory{k} stores the reciprocal of the inner product between
% sHistory{k} and yHistory{k}.
rhoHistory = cell(1, options.memory);

% Scaling of direction given by getDirection for acceptable step
alpha = 1;

% Scaling of initial matrix, Barzilai-Borwein.
scaleFactor = 1;

% Norm of the step
stepsize = 1;

% Stores whether the step is accepted by the cautious update check.
accepted = true;

% Compute cost function and its gradient
[hCurCost,hCurGradient] = alignment_costgrad(q1,q2tilde,M);
hCurGradNorm = M.norm(hCurGradient);

% Line-search statistics for recording in info.
lsstats = [];

% Flag to control restarting scheme to avoid infinite loops (see below)
ultimatum = false;

% Save stats in a struct array info, and preallocate.
stats = savestats();
info(1) = stats;
info(min(10000,options.maxiter+1)).iter = [];

if options.verbosity >= 2
    fprintf(' iter                   cost val            grad. norm           alpha\n');
end

% Plot the evolution of q2 being iteratively warped to optimally
% match q1.
if options.plotevol
    figure(100); clf; plot(M.t,q1,'b',M.t,q2tilde,'r');
    %         figure(101); clf; plot(M.t,htilde);
    pause;
end

% Main iteration
while true
    
    % Display iteration information
    if options.verbosity >= 2
        fprintf('%5d    %+.16e        %.8e      %.4e\n', ...
            k, hCurCost, hCurGradNorm, alpha);
    end
    
    % Start timing this iteration
    timetic = tic();
    
    % Run standard stopping criterion checks
    [stop, reason] = stoppingcriterion(options, info, k+1);
    
    % If none triggered, run specific stopping criterion check
    if ~stop
        if stats.stepsize < options.minstepsize
            % To avoid infinite loop and to push the search further
            % in case BFGS approximation of Hessian is off towards
            % the end, we erase the memory by setting k = 0;
            % In this way, it starts off like a steepest descent.
            % If even steepest descent does not work, then it is
            % hopeless and we will terminate.
            if ~ultimatum
                if options.verbosity >= 2
                    fprintf(['stepsize is too small, restarting ' ...
                        'the bfgs procedure at the current point.\n']);
                end
                j = 0;
                ultimatum = true;
            else
                stop = true;
                reason = sprintf(['Last stepsize smaller than '  ...
                    'minimum allowed; options.minstepsize = %g.'], ...
                    options.minstepsize);
            end
        else
            % We are not in trouble: lift the ultimatum if it was on.
            ultimatum = false;
        end
    end
    
    if stop
        if options.verbosity >= 1
            fprintf([reason '\n']);
        end
        break;
    end
    
    % Compute BFGS direction
    p = getDirection(M,hCurGradient,sHistory,yHistory,rhoHistory,...
        scaleFactor,min(j,options.memory));
    
    % Execute line-search
    [stepsize,hNext,lsstats] = linesearch_hint(M,p,hCurCost,...
        M.inner(hCurGradient,p),q1,q2tilde,options);
    
    % Iterative update of optimal diffeomorphism and q2 via group action
    htilde=group_action_SRVF(htilde,hNext,M);
    q2tilde=group_action_SRVF(q2tilde,hNext,M);
    
    % Record the BFGS step-multiplier alpha which was effectively
    % selected. Toward convergence, we hope to see alpha = 1.
    alpha = stepsize/M.norm(p);
    step = alpha*p;
    
    % Query cost and gradient at the candidate new point.
    [hNextCost,hNextGradient] = alignment_costgrad(q1,q2tilde,M);
    
    % Compute sk and yk
    sk = step;
    yk = hNextGradient-hCurGradient;
    
    % Computation of the BFGS step is invariant under scaling of sk and
    % yk by a common factor. For numerical reasons, we scale sk and yk
    % so that sk is a unit norm vector.
    norm_sk = M.norm(sk);
    sk = sk/norm_sk;
    yk = yk/norm_sk;
    
    inner_sk_yk = M.inner(sk,yk);
    inner_sk_sk = M.norm(sk)^2;    % ensures nonnegativity
    
    
    % If the cautious step is accepted (which is the intended
    % behavior), we record sk, yk, and rhok and need to do some
    % housekeeping. If the cautious step is rejected, these are not
    % recorded. In all cases, hNext is the next iterate: the notion of
    % accept/reject here is limited to whether or not we keep track of
    % sk, yk, rhok to update the BFGS operator.
    cap = options.strict_inc_func(hCurGradNorm);
    if inner_sk_sk~=0 && (inner_sk_yk/inner_sk_sk)>=cap
        
        accepted = true;
        
        rhok = 1/inner_sk_yk;
        
        scaleFactor = inner_sk_yk/M.norm(yk)^2;
        
        % Time to store the vectors sk, yk and the scalar rhok.
        
        % If we are out of memory
        if j>=options.memory
            
            % sk and yk are saved from 1 to the end with the most
            % current recorded to the rightmost hand side of the cells
            % that are occupied. When memory is full, do a shift so
            % that the rightmost is earliest and replace it with the
            % most recent sk, yk.
            if options.memory>1
                sHistory = sHistory([2:end,1]);
                yHistory = yHistory([2:end,1]);
                rhoHistory = rhoHistory([2:end,1]);
            end
            if options.memory>0
                sHistory{options.memory} = sk;
                yHistory{options.memory} = yk;
                rhoHistory{options.memory} = rhok;
            end
            
            % If we are not out of memory
        else
            
            sHistory{j+1} = sk;
            yHistory{j+1} = yk;
            rhoHistory{j+1} = rhok;
            
        end
        
        j = j+1;
        
        % The cautious step is rejected: we do not store sk, yk, rhok.
    else
        
        accepted = false;
        
    end
    
    % Update variables to new iterate.
    hCurGradient = hNextGradient;
    hCurGradNorm = M.norm(hNextGradient);
    hCurCost = hNextCost;
    
    % Plot the evolution of q2 being iteratively warped to optimally
    % match q1.
    if options.plotevol
        figure(100); clf; plot(M.t,q1,'b',M.t,q2tilde,'r');
        %         figure(101); clf; plot(M.t,hk);
        pause;
    end
    
    % iter is the number of iterations we have accomplished.
    k = k+1;
    
    % Log statistics for freshly executed iteration
    stats = savestats();
    info(k+1) = stats;
    
end

% Housekeeping before we return
info = info(1:k+1);
gammaOpt = cumtrapz(M.t,htilde.^2);
gammaOpt = gammaOpt/gammaOpt(end);
q2Opt = q2tilde;
cost = hCurCost;

if options.verbosity >= 1
    fprintf('Total time is %f [s] (excludes statsfun)\n', ...
        info(end).time);
end


% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = k;
        stats.cost = hCurCost;
        stats.gradnorm = hCurGradNorm;
        if k==0
            stats.stepsize = NaN;
            stats.time = toc(timetic);
            stats.accepted = NaN;
        else
            stats.stepsize = stepsize;
            stats.time = info(k).time + toc(timetic);
            stats.accepted = accepted;
        end
        stats.linesearch = lsstats;
        %         stats = applyStatsfun(problem, xCur, storedb, key, options, stats);
    end

end


% BFGS step, see Wen's paper for details. This functon takes in a tangent
% vector g, and applies an approximate inverse Hessian P to it to get Pg.
% Then, -Pg is returned. Parallel transport is not needed for this problem
% since we always work in the tangent space of identity.
function dir = getDirection(M,hCurGradient,sHistory,yHistory,rhoHistory,...
    scaleFactor,j)

q=hCurGradient;
inner_s_q=zeros(1,j);

for i=j:-1:1
    inner_s_q(1,i)=rhoHistory{i}*M.inner(sHistory{i},q);
    q=q-inner_s_q(1,i)*yHistory{i};
end

r=scaleFactor*q;

for i=1:j
    omega=rhoHistory{i}*M.inner(yHistory{i},r);
    r=r+(inner_s_q(1,i)-omega)*sHistory{i};
end

dir=-r;

end

function opts = mergeOptions(opts_sub, opts_master)
% Merges two options structures with one having precedence over the other.
%
% function opts = mergeOptions(opts1, opts2)
%
% input: opts1 and opts2 are two structures.
% output: opts is a structure containing all fields of opts1 and opts2.
% Whenever a field is present in both opts1 and opts2, it is the value in
% opts2 that is kept.
%
% The typical usage is to have opts1 contain default options and opts2
% contain user-specified options that overwrite the defaults.
%
% See also: getGlobalDefaults

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors:
% Change log:


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

function f = alignment_cost(h,q1,q2k,M)

% Evaluate the cost function f = ||q1 - ((q2,hk),h)||^2.
% h=sqrt{\dot{\gamma}} is a sequential update of cumulative warping hk

t=M.t;
q2new=group_action_SRVF(q2k,h,M);
f=normL2(q1-q2new,t)^2; % Cost

end

function [f,g] = alignment_costgrad(q1,q2k,M)

% Evaluate the cost function f = ||q1 - (q2,hk)||^2, and
% evaluate the gradient g = grad f in the tangent space of identity.
% hk=sqrt{\dot{\gamma_k}} is the cumulative warping of q2 produced by an
% iterative sequential optimization algorithm.

t=M.t;
T=M.T;

% Compute cost
f=normL2(q1-q2k,t)^2;

% Compute cost gradient
q2kdot=gradient(q2k,1/(T-1));
dq=q1-q2k;
v=2*cumtrapz(t,sum(dq.*q2kdot,1))-sum(dq.*q2k,1);
g=v-trapz(t,v);
end

function val = normL2(f,t)

val=sqrt(innerProdL2(f,f,t));
end

function [stop, reason] = stoppingcriterion(options, info, last)
% Checks for standard stopping criteria, as a helper to solvers.
%
% function [stop, reason] = stoppingcriterion(problem, x, options, info, last)
%
% Executes standard stopping criterion checks, based on what is defined in
% the info(last) stats structure and in the options structure.
%
% The returned number 'stop' is 0 if none of the stopping criteria
% triggered, and a (strictly) positive integer otherwise. The integer
% identifies which criterion triggered:
%  0 : Nothing triggered;
%  1 : Cost tolerance reached;
%  2 : Gradient norm tolerance reached;
%  3 : Max time exceeded;
%  4 : Max iteration count reached;
%  6 : User defined stopfun criterion triggered.
%
% The output 'reason' is a string describing the triggered event.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors:
% Change log:
%
%   Apr. 2, 2015 (NB):
%       'reason' now contains the option (name and value) that triggered.
%
%   Aug. 3, 2018 (NB):
%       Removed check for costevals, as it was never used, and the new
%       manopt counters allow to do this in a more transparent way.
%       Furthermore, now, options.stopfun can have 1 or 2 outputs: the
%       first is a boolean indicating whether or not to stop, and the
%       (optional) second output is a string indicating the reason.


stop = 0;
reason = '';

stats = info(last);

% Target cost attained
if isfield(stats, 'cost') && isfield(options, 'tolcost') && ...
        stats.cost <= options.tolcost
    reason = sprintf('Cost tolerance reached; options.tolcost = %g.', options.tolcost);
    stop = 1;
    return;
end

% Target gradient norm attained
if isfield(stats, 'gradnorm') && isfield(options, 'tolgradnorm') && ...
        stats.gradnorm < options.tolgradnorm
    reason = sprintf('Gradient norm tolerance reached; options.tolgradnorm = %g.', options.tolgradnorm);
    stop = 2;
    return;
end

% Allotted time exceeded
if isfield(stats, 'time') && isfield(options, 'maxtime') && ...
        stats.time >= options.maxtime
    reason = sprintf('Max time exceeded; options.maxtime = %g.', options.maxtime);
    stop = 3;
    return;
end

% Allotted iteration count exceeded
if isfield(stats, 'iter') && isfield(options, 'maxiter') && ...
        stats.iter >= options.maxiter
    reason = sprintf('Max iteration count reached; options.maxiter = %g.', options.maxiter);
    stop = 4;
    return;
end

%     % Check whether the possibly user defined stopping criterion
%     % triggers or not.
%     if isfield(options, 'stopfun')
%         % options.stopfun can have 1 or 2 outputs, but checking this with
%         % nargout does not always work because it is technical to determine
%         % for anonymous functions. Thus, we use our best guess. Nargout
%         % returns -1 when it cannot determine the number of outputs, in
%         % which case we take the safer approach of assuming 1 output.
%         switch nargout(options.stopfun)
%             case 2
%                 [userstop, reason] = options.stopfun(problem, x, info, last);
%             case {1, -1}
%                 userstop = options.stopfun(problem, x, info, last);
%                 reason = ['User defined stopfun criterion triggered; ' ...
%                           'see options.stopfun.'];
%             otherwise
%                 error('manopt:stoppingcriterion:stopfunoutputs', ...
%                       'options.stopfun must have one or two outputs.');
%         end
%         if userstop
%             stop = 6;
%             if nargout(options.stopfun) == -1
%                 reason = [reason, '\n(A reason may have been ' ...
%                           'provided, but stoppingcriterion was ' ...
%                           'unable to determine\nthe number of ' ...
%                           'output arguments of options.stopfun.)'];
%             end
%             return;
%         end
%     end

end

function [stepsize,newh,lsstats] = linesearch_hint(M,d,f0,df0,q1,q2k,options)

% Armijo line-search based on the line-search hint in the problem structure.
%
% Base line-search algorithm for descent methods, based on a simple
% backtracking method. The search direction provided has to be a descent
% direction, as indicated by a negative df0 = directional derivative of f
% at the identity element along d.
%
% The algorithm selects a hardcoded initial step size. If that
% step does not fulfill the Armijo sufficient decrease criterion, that step
% size is reduced geometrically until a satisfactory step size is obtained
% or until a failure criterion triggers.
%
% Below, the step is constructed as alpha*d, and the step size is the norm
% of that vector, thus: stepsize = alpha*norm_d. The step is executed by
% computing the exponential mapping exp_{hid}(alpha*d), giving newh.
%
% Inputs/Outputs : see help for linesearch
%
% See also: steepestdescent conjugategradients linesearch

% This file is a modified version of linesearch_hint.m in Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 17, 2014.
% Modified by: Darshan Bryner, August 1, 2019.


% Backtracking default parameters. These can be overwritten in the
% options structure which is passed to the solver.
default_options.ls_contraction_factor = .5;
default_options.ls_suff_decr = 1e-6;
default_options.ls_max_steps = 25;
default_options.ls_backtrack = true;
default_options.ls_force_decrease = true;

if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(default_options, options);

contraction_factor = options.ls_contraction_factor;
suff_decr = options.ls_suff_decr;
max_ls_steps = options.ls_max_steps;

% Initialize alpha.
alpha=1;

% Identity element, i.e. where we are located on the manifold.
hid=ones(1,M.T);

% Make the chosen step and compute the cost there.
newh = M.exp(hid, d, alpha);
newf = alignment_cost(newh,q1,q2k,M);
cost_evaluations = 1;

% Backtrack while the Armijo criterion is not satisfied or if newh goes outside positive orthant.
while options.ls_backtrack && (newf > f0 + suff_decr*alpha*df0 || sum(newh<=0)>0)
    
    % Reduce the step size,
    alpha = contraction_factor * alpha;
    
    % and look closer down the line.
    newh = M.exp(hid, d, alpha);
    newf = alignment_cost(newh,q1,q2k,M);
    cost_evaluations = cost_evaluations + 1;
    
    % Make sure we don't run out of budget.
    if cost_evaluations >= max_ls_steps
        break;
    end
    
end

% If we got here without obtaining a decrease, we reject the step.
if options.ls_force_decrease && newf > f0
    alpha = 0;
    newh = hid;
    newf = f0; %#ok<NASGU>
end

% As seen outside this function, stepsize is the size of the vector we
% retract to make the step from h to newh. Since the step is alpha*d:
norm_d = M.norm(d);
stepsize = alpha * norm_d;

% Return some statistics also, for possible analysis.
lsstats.costevals = cost_evaluations;
lsstats.stepsize = stepsize;
lsstats.alpha = alpha;

end

function [qnew,gamma] = group_action_SRVF(q,h,M)

    % Computes the diffeomorphis group action on an SRVF given by 
    % (q,\gamma) = q(\gamma(t))h(t), where h=\sqrt{\dot{\gamma}}. 
    
    [p,~]=size(q);
    gamma=cumtrapz(M.t,h.^2);
    gamma=gamma/gamma(end);
    h=sqrt(gradient(gamma,M.t));
    if p>1
        h=repmat(h,p,1);
    end
    qnew=spline(M.t,q,gamma).*h;
end