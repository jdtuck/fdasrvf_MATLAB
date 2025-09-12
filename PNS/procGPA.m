function [rotated, ProcrustesMean, centroidsize] = procGPA(data,reflection)
% procGPA : PROCedure for full Generalized Procrustes Analysis
%           using M function 'procrustes'
%
%   [rotated, ProcrustesMean, centroidsize] = procGPA(data)
%                   with k x m x n matrix of landmark data.
%                     (reflection is prohibited)
%   [rotated, ProcrustesMean, centroidsize] = procGPA(data,reflection)
%                   with k x m x n matrix of landmark data,
%                       reflection = true or false
%                       to indicate if reflection is included
%
%    k : number of landmarks
%    m : dimension of space that landmarks lie
%    n : number of samples
%
% Jan 15, 2010
% Sungkyu Jung


        
ireflection = false;
if nargin ==2
    ireflection = reflection;
end
[k m n] =size(data);
H = Helmertsub(k);
C = H'*H; % centering matrix

% step 1. centering and obtaining the centroidsize of each configuration
centroidsize = zeros(n,1);
for i = 1:n;
    centeredconfigi = C*data(:,:,i);
    centroidsize(i) = norm(centeredconfigi,'fro');
    data(:,:,i) = centeredconfigi;
end



% step 1.5 initial rotation including scaling;
for i = 1:n;
    xbari = sum(data(:,:,[1:(i-1) (i+1):end]),3)/(n-1);
    if ireflection
        [d, xi] = procrustes(xbari, data(:,:,i),'reflection','best','scaling',true);
    else
        [d, xi] = procrustes(xbari, data(:,:,i),'reflection',false,'scaling',true);
    end
    data(:,:,i) = xi;
end
diff1 = 1;
Gnow = ssnpd(data);
while diff1 > 1e-15
    
    % step 2. For each ith configuration, superimpose xi to the xbari repeat
    % until there is no change;
    % criteria: ssnpd - sum of squared norms of pairwise difference
    
    diff = 1;
    while diff > 1e-15;
        for i = 1:n;
            xbari = sum(data(:,:,[1:(i-1) (i+1):end]),3)/(n-1);
            if ireflection
                [d, xi] = procrustes(xbari, data(:,:,i),'reflection','best','scaling',false);
            else
                [d, xi] = procrustes(xbari, data(:,:,i),'reflection',false,'scaling',false);
            end
            data(:,:,i) = xi;
        end
        Gnext = ssnpd(data);
        diff = Gnow - Gnext;
        Gnow = Gnext;
    end
    % step 3. scaling
    
    vecdata = reshape(data,k*m,n);
    RHO = corr(vecdata);
    [V,D] = eig(RHO);
    phi = V(:,end);
    xipnorm2arr = sum(vecdata.^2);
    sumxipnorm2arr = sum(xipnorm2arr);
    for i=1:n
        betai = phi(i)*sqrt(sumxipnorm2arr/xipnorm2arr(i));
        data(:,:,i) = betai*data(:,:,i);
    end
    Gnext = ssnpd(data);
    diff1 = Gnow - Gnext;
    Gnow = Gnext;
end
rotated = data;
ProcrustesMean = mean(data,3);
return  ;


function g = ssnpd(data)
[k d n] = size(data);
g = 0;
for i = 1:n;
    for j = (i+1):n
        xi = data(:,:,i);
        xj = data(:,:,j);
        g = g + norm(xi-xj,'fro')^2;
    end
end
g = g / n ;
return
