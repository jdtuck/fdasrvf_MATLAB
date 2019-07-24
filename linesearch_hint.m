function [stepsize, newh, lsstats] = linesearch_hint(M, d, f0, iter, df0, q1, q2, lambda, h0, options)
    % Armijo line-search based on the line-search hint in the problem structure.
    %
    % function [stepsize, newh, lsstats] = 
    %            linesearch_hint(M, d, f0, iter, df0, q1, q2, options)
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
        
        % Initialize alpha. For elastic alignment, it seems the first iteration
        % can be a bit unstable, so start with alpha=0.01. 
        if iter==0
            alpha=0.01;
        else
            alpha=1;
        end
        
        % Identity element, i.e. where we are located on the manifold.
        hid=ones(1,M.T);
    
        % Make the chosen step and compute the cost there.
        newh = M.exp(hid, d, alpha);
        newf = alignment_cost(newh,h0,q1,q2,M,lambda);
        cost_evaluations = 1;
        
        % Backtrack while the Armijo criterion is not satisfied.
        while options.ls_backtrack && newf > f0 + suff_decr*alpha*df0
            
            % Reduce the step size,
            alpha = contraction_factor * alpha;
            
            % and look closer down the line.
            newh = M.exp(hid, d, alpha);
            newf = alignment_cost(newh,h0,q1,q2,M,lambda);
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