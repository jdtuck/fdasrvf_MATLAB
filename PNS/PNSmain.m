function [resmat PNS]= PNSmain(data,itype,alpha,R)
% PNSMAIN Analysis of Principal Nested Spheres for data on hyperspheres
% [resmat PNS]= PNSmain(data)
% [resmat PNS]= PNSmain(data,itype)
% [resmat PNS]= PNSmain(data,itype,alpha)
% [resmat PNS]= PNSmain(data,itype,alpha,R)
%
% Input:
%   data     (d+1) x n data matrix where each column is a unit vector.
%
%   itype    0  'seq.test' : (default) ordinary Principal Nested Sphere
%                               with sequential tests.
%            1  'small'    : Principal Nested SMALL Sphere
%            2  'great'    : Principal Nested GREAT Sphere (radius pi/2)
%
%  alpha     0.05 (default) : size of Type I error allowed for each test,
%            could be any number between 0 and 1.
%  R         100 (default) : number of bootsrap samples to be evaluated for
%            the sequential test.
%
% Output:
%   resmat  The commensurate residual matrix (X_PNS). Each entry in row k
%           works like the kth principal component score.
%
%   PNS     Matlab structure array with following fields:
%       mean     : location of the PNSmean
%       radii    : size (radius) of PNS
%       orthoaxis: orthogonal axis 'v_i' of subspheres
%       dist     : distance 'r_i' of subspheres
%       pvalues  : p-values of LRT and parametric boostrap tests (if any)
%       ratio    : estimated ratios
%       itype    : type of methods for fitting subspheres


if nargin == 1;
    itype = 0; %  'seq.test'; % default: sequential test
    alpha = 0.05;
    R = 100;
end

if nargin == 2;
    alpha = 0.05;
    R = 100;
end

if nargin == 3;
    R = 100;
end


[k n] = size(data); % data on (k-1)-sphere in Real^k
[uu ll]=svd(data);
maxd = find(diag(ll) < 1e-15,1);
if isempty(maxd) || k > n
    maxd = min(k,n)+1;
end
nullspdim = k - maxd + 1; % dimension of subspace that contains no data

d = k - 1; % dimension of the sphere
disp(['Message from PNSmain.m; dataset is on ' num2str(d) '-sphere.']);

if nullspdim > 0
    disp([' .. found null space of dimension ' num2str(nullspdim)...
        ',to be trivially reduced.']);
end
resmat = zeros(d,n); % d dimensional residual matrix
% there will be d-1 subspheres
orthaxis = cell(d-1,1);
dist = zeros(d-1,1);
pvalues = zeros(d-1,2);
ratio = zeros(d-1,1);

% (HDLSS case) fit nested great spheres for dimension reduction
% where no residual is present.
currentSphere = data;
for i = 1:nullspdim;
    oaxis = uu(:,end-i+1); % orthogonal axis for subsphere
    r = pi/2;              % distance for subsphere
    pvalues(i,:) = [NaN,NaN];        % No test performed (p-value is Not-a-Number)
    res = acos(oaxis'*currentSphere)-r; % residuals
    % save subsphere parameters
    orthaxis{i} = oaxis; dist(i) = r;
    % save residuals
    resmat(i,:) = res;
    % projection to subsphere and transformation to isomorphic sphere
    NestedSphere = rotMat(oaxis)*currentSphere;
    currentSphere = NestedSphere(1:(k-i),:)./...
        repmat(sqrt(1-NestedSphere(end,:).^2),k-i,1);
    % transform singular vectors accordingly
    uu = rotMat(oaxis)*uu;
    uu = uu(1:(k-i),:)./repmat(sqrt(1-uu(end,:).^2),k-i,1); 
end 

if itype == 0;  % 'seq.test'
    disp([' .. sequential tests with significance level ' num2str(alpha)])
    isIsotropic = false; % to be used as a test procedure
    for i=nullspdim+1:(d-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0); % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs;
                [centerg, rg] = getSubSphere(currentSphere,1); % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg;
                
                % Chi-squared statistic for a likelihood test
                pval1 = LRTpval(resGREAT,resSMALL,n);
                
                pvalues(i,1) = pval1;
                if pval1 > alpha
                    center = centerg; r = rg;
                    pvalues(i,2) = NaN;
                    disp([num2str(d-i+1) '-sphere to '...
                        num2str(d-i) '-sphere, by '...
                        'GREAT sphere' ', p(LRT) = ' num2str(pval1)]);
                else
                    pval2 = vMFtest(currentSphere,R);
                    pvalues(i,2) = pval2;
                    if pval2 > alpha
                        center = centerg; r = rg;
                        disp([num2str(d-i+1) '-sphere to '...
                        num2str(d-i) '-sphere, by '...
                            'GREAT sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(vMF) = ' num2str(pval2)]);
                        isIsotropic = true;
                    else
                        center = centers; r = rs;
                        disp([num2str(d-i+1) '-sphere to '...
                        num2str(d-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(vMF) = ' num2str(pval2)]);            
                    end
                end
                    
            case true
                [center, r] = getSubSphere(currentSphere,1);
                disp([num2str(d-i+1) '-sphere to '...
                        num2str(d-i) '-sphere, by '...
                    'GREAT sphere, restricted by testing vMF distn']);
                pvalues(i,1) = NaN;
                pvalues(i,2) = NaN;
        end
        res = acos(center'*currentSphere)-r;
        % save subsphere parameters
        orthaxis{i} = center; dist(i) = r;
        % save residuals
        resmat(i,:) = res;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere;
        currentSphere = NestedSphere(1:(k-i),:)./...
            repmat(sqrt(1-NestedSphere(end,:).^2),k-i,1);
    end
elseif itype == 1|| itype ==2 % strictly 'small' or 'great' spheres
    pvalues = NaN;
    for i=nullspdim+1:(d-1)
        % estimate the best fitting subsphere
        %  with small sphere if itype = 1
        %  with  great sphere if itype = 2
        [center, r] = getSubSphere(currentSphere,itype-1);
        res = acos(center'*currentSphere)-r;
        % save subsphere parameters
        orthaxis{i} = center; dist(i) = r;
        % save residuals
        resmat(i,:) = res;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere;
        currentSphere = NestedSphere(1:(k-i),:)./...
            repmat(sqrt(1-NestedSphere(end,:).^2),k-i,1);
    end
else
    disp('!!! Error from PNSmain.m:');
    disp('!!! itype must be 0 (seq.test), 1 (small), or 2 (great)');
    disp('!!!   Terminating execution     ') ;
    return ;
end

%
% scatter(currentSphere(1,:),currentSphere(2,:)); axis equal;
% hold on; scatter(cos(meantheta),sin(meantheta),'r')
% scatter(cos(S1toRadian),sin(S1toRadian),'.k'); axis equal;

% currentSphere has (intrinsic) dimension 1
% compute PNSmean and deviations.

% parametrize 1-sphere to angles
S1toRadian = atan2(currentSphere(2,:),currentSphere(1,:));
% Geodesic mean of angles
meantheta = geodmeanS1(S1toRadian');
orthaxis{d} = meantheta;
% save deviations from PNSmean
resmat(d,:) =mod(S1toRadian - meantheta + pi,2*pi)-pi;


radii =1;
for i = 1:(d-1)
    radii = [radii; prod(sin(dist(1:i)))];
end
resmat = flipud(repmat(radii,1,n).*resmat); % scaled residuals

PNS.radii = radii;         % size (radius) of nested spheres from largest to smallest

PNS.orthaxis = orthaxis;    % orthogonal axis of (d-1) subspheres and the anglemean for PNSmean
PNS.dist = dist;            % distances for (d-1) subspheres
PNS.pvalues = pvalues;      % d-1 pvalues from sequential tests
PNS.ratio = ratio;   % d-1 ratios estimated
PNS.basisu = [];
PNS.mean = PNSe2s(zeros(d,1),PNS); % PNSmean of the data

switch itype
    case 0
        PNS.itype = 'seq.test';
    case 1
        PNS.itype = 'small';
    case 2
        PNS.itype = 'great';
end

end

function pval = LRTpval(resGREAT,resSMALL,n)
chi2 = max(n*log(sum(resGREAT.^2)/sum(resSMALL.^2)),0);
pval = 1-chi2cdf(chi2,1); % p-value of the likelihood test
end