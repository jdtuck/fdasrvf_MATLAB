function q = f_to_srvf(f,time,smooth,spl, parallel)
% F_TO_SRVF Convert function to Square-Root Velocity Function
% -------------------------------------------------------------------------
% Convert to SRSF
%
% Usage: q = f_to_srvf(f,time)
%
% This function converts functions to srsf
%
% Input:
% f: matrix of functions
% time: vector of time samples
% smooth: use smoothing splines (default: true)
% spl: use b-spline computation (old behavior)
% paralell: compute in parallel
%
% Output:
% q: matrix of SRSFs

arguments
    f double
    time double
    smooth=true
    spl=false
    parallel=false
end

if parallel == 1
    if isempty(gcp('nocreate'))
        % prompt user for number threads to use
        nThreads = input('Enter number of threads to use: ');
        if nThreads > 1
            parpool(nThreads);
        elseif nThreads > 12 % check if the maximum allowable number of threads is exceeded
            while (nThreads > 12) % wait until user figures it out
                fprintf('Maximum number of threads allowed is 12\n Enter a number between 1 and 12\n');
                nThreads = input('Enter number of threads to use: ');
            end
            if nThreads > 1
                parpool(nThreads);
            end
        end
    end
end

binsize = mean(diff(time));
[M, N] = size(f);

if smooth
    fy = zeros(M,N);
    if parallel
        parfor ii = 1:N
            y = fit(time(:), f(:,ii),'smoothingspline');
            fy(:,ii) = differentiate(y, time);
        end
    else
        for ii = 1:N
            y = fit(time(:), f(:,ii),'smoothingspline');
            fy(:,ii) = differentiate(y, time);
        end
    end
    q = fy./sqrt(abs(fy)+eps);
elseif spl
    fy = zeros(M,N);
    if parallel
        parfor ii = 1:N
            y = Bspline(f(:,ii),3);
            ydiff = diff(y);
            fy(:,ii) = ydiff(1:length(time))/binsize;
            f(:,ii) = y(1:length(time));
        end
    else
        for ii = 1:N
            y = Bspline(f(:,ii),3);
            ydiff = diff(y);
            fy(:,ii) = ydiff(1:length(time))/binsize;
            f(:,ii) = y(1:length(time));
        end
    end
    q = fy./sqrt(abs(fy)+eps);
else
    if size(f,2)>1
        [~,fy] = gradient(f,binsize,binsize);
    else
        fy = gradient(f,binsize);
    end
    q = fy./sqrt(abs(fy)+eps);
end
