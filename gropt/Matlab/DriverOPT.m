function [output, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, InitialX)
    for i = 1 : length(ManiParams)
        if(~isfield(ManiParams(i), 'numofmani'))
            ManiParams(i).numofmani = 1;
        end
        if(~isfield(ManiParams(i), 'n') || isempty(ManiParams(i).n))
            ManiParams(i).n = 1;
        end
        if(~isfield(ManiParams(i), 'm') || isempty(ManiParams(i).m))
            ManiParams(i).m = 1;
        end
        if(~isfield(ManiParams(i), 'p') || isempty(ManiParams(i).p))
            ManiParams(i).p = 1;
        end
        if(~isfield(ManiParams(i), 'ParamSet') || isempty(ManiParams(i).ParamSet))
            ManiParams(i).ParamSet = 1;
        end
    end
    
    global xsizeparams
    xsizeparams.main = size(InitialX.main);
    fh = @(x)f(x, fhandle);
    Egfh = @(x)Egf(x, gfhandle);
    EHh = @(x, eta)EHess(x, eta, Hesshandle);

    [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = DriverMexProb(fh, Egfh, EHh, SolverParams, ManiParams, HasHHR, InitialX);
    output.main = reshape(FinalX.main, xsizeparams.main);
end

function [output, x] = f(x, fhandle)
    global xsizeparams
    x.main = reshape(x.main, xsizeparams.main);
    [output, x] = fhandle(x);
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), [], 1);
    end
end

function [output, x] = Egf(x, gfhandle)
    global xsizeparams
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    [output, x] = gfhandle(x);
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), [], 1);
    end
end

function [output, x] = EHess(x, eta, Hesshandle)
    global xsizeparams
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    % eta.main has the same size as x.main
    eta.main = reshape(eta.main, xsizeparams.main);
    [output, x] = Hesshandle(x, eta);
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), [], 1);
    end
end
