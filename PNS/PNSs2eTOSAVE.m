function EuclidData= PNSs2e(spheredata,PNS,ialign)
% PNSS2E : PNS Sphere to Euclidean-type representation
% Coordinates in Euclidean type representation by PNS from points in sphere
%
%   EuclidData = PNSs2e(A,PNS) with d x n matrix A consists of column
%       vectors that are on the (d-1) sphere, and PNS is the output from
%       PNSmain.
%   EuclidData = PNSs2e(spheredata,PNS,1) to force align spheredata to
%                                       PNS.alignbase
%

if size(spheredata,1)~=size(PNS.mean,1)
    disp(' Error from PNSs2e.m ')
    disp(' Dimensions of the sphere and PNS decomposition do not match');
    return;
end

[kk n] = size(spheredata);

iAlignFlag = false;
if nargin > 2;
    iAlignFlag = (ialign == 1);
end
if iAlignFlag % align to PNS.alignbase
    [k m]=size(PNS.alignbase);
    for i=1:n
        rotZ = optrotPreshape(reshape(spheredata(:,i),k,m),PNS.alignbase);
        spheredata(:,i) = rotZ(:);
    end
end

Res = zeros(kk-1,n); % k-1 dimensional residual matrix
currentSphere = spheredata;

for i = 1:(kk-2);
    v = PNS.orthaxis{i};
    r = PNS.dist(i);
    res = acos(v'*currentSphere)-r;
    Res(i,:) = res;
    NestedSphere = rotMat(v)*currentSphere;
    currentSphere = NestedSphere(1:(kk-i),:)./...
        repmat(sqrt(1-NestedSphere(end,:).^2),kk-i,1);
end

S1toRadian = atan2(currentSphere(2,:),currentSphere(1,:));
devS1 = mod(S1toRadian - PNS.orthaxis{end} + pi,2*pi)-pi;
Res(kk-1,:) =devS1;

EuclidData = flipud(repmat(PNS.radii,1,n).*Res); % scaled residuals
