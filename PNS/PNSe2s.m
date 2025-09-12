function T= PNSe2s(resmat,PNS)
% PNSE2S PNS coordinate transform from Euclidean-type residual matrix to Sphere
% Spheredata = PNSe2s(data,PNS)
%      where 'data' is d x m data matrix in PNS coordinate system (for any
%      m >= 1), 'PNS' is the structural array from PNSmain.m
%
%

[dm, n] = size(resmat);
% dm is the intrinsic dimension of the sphere
%    or the reduced sphere in HDLSS case.
% n  is the sample size.
NSOrthaxis = flipud(PNS.orthaxis(1:end-1));
NSradius = flipud(PNS.dist);
geodmean = PNS.orthaxis{end};

% "standardize" the coordinates
res = resmat./repmat(flipud(PNS.radii),1,n);

% iteratively mapping back to S^d

% S^1 to S^2
T = rotMat(NSOrthaxis{1})'* ...
    [repmat(sin(NSradius(1)+res(2,:)),2,1).*[cos(geodmean+res(1,:)) ; sin(geodmean+res(1,:))] ;
    cos(NSradius(1)+res(2,:))];
% S^2 to S^d
for i=1:dm-2
    T = rotMat(NSOrthaxis{i+1})'*...
        [repmat(sin(NSradius(i+1)+res(i+2,:)),2+i,1).*T ;
        cos(NSradius(i+1)+res(i+2,:))];
end

if isempty(PNS.basisu) == 0
    % Then this is the HDLSS case
    T = PNS.basisu*T;
end

% check the entries of the following are all one
% sum(T.^2)