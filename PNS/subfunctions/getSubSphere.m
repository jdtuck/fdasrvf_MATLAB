function [center, r] = getSubSphere(data,geodesic)
% GETSUBSPHERE : The least square estimates of the best fitting subsphere 
%                to the data on the unit hyper-sphere. 
% [center]= getSubSphere(data), with d x n data matrix with 
%           each column having unit length, returns the center
% [center, r]= getSubSphere(data), with d x n data matrix with each 
%              column having unit length, returns the center and the
%              geodesic radius.
% [center, r]= getSubSphere(data,1) forces the subsphere radius as 1
%
% % example
% n =50;
% theta = linspace(0,pi*1.5,n);
% data = ExpNP([cos(theta); sin(theta)] + 1*randn(2,n));
% data = [data ; zeros(1,n)];
% data = data + randn(size(data));
% [d n ] = size(data);
% data = data./repmat(sqrt(sum(data.^2)),d, 1);
% 
%    See also NestedSpheres, LMFSPHEREFit.

% Last updated Sep 20, 2009
% Sungkyu Jung

if nargin==0 && nargout==0, help getSubSphere, return, end     %   Display help

if nargin==1
   geodesic = 0;  % Do not force the radius to be pi/2
end

% 
% % % Choice of initial values
% % last singular vector
% initialCenter = initilalc(data);
% 
% c0 = initialCenter;
% TOL = 1e-10;
% cnt =0;
% err = 1;
% [d n]= size(data);
% Gnow = 1e+10;
% 
% while err > TOL
%     c0 = c0/norm(c0); % normalize the new candidate
%     rot = rotMat(c0);  % rotation matrix : c0 -> North Pole
%     TpData = LogNPd(rot*data); % Tangent projection by Log map 
%     [newCenterTp r]= LMFsphereFit(TpData,zeros(d-1,1),geodesic);  % Fit a circle at Tp
%     if r > pi
%         disp(['Message from getSubSphere.m: (circfit) large radius = ' num2str(r) ]);
%         r = pi/2;
%         [U dd] = svd(TpData);
%         newCenterTp = U(:,end)*pi/2;
%     end
%     newCenter = ExpNPd(newCenterTp); % Bring back to the sphere by Exp map
%     center = rot\newCenter;  % (inv(rot)*newCenter) rotate back the newCenter
%     Gnext = objfn(center,r,data);
%     err = abs(Gnow - Gnext) ;
%     Gnow = Gnext;
%     c0 = center;
%     cnt = cnt+1;
%     if cnt > 30;
%        disp(['Message from getSubSphere.m: (circfit) iteration reached 30th with error ' num2str(err)]);
%        break;
%     end
% end
% i2save.Gnow = Gnow;
% i2save.center= center;
% i2save.r = r;


% last singular vector
[U, dd] = svd(data);
initialCenter = U(:,end);

c0 = initialCenter;
TOL = 1e-10;
cnt =0;
err = 1;
[d, n]= size(data);
Gnow = 1e+10;

while err > TOL
    c0 = c0/norm(c0); % normalize the new candidate
    rot = rotMat(c0);  % rotation matrix : c0 -> North Pole
    TpData = LogNPd(rot*data); % Tangent projection by Log map 
    [newCenterTp, r]= LMFsphereFit(TpData,zeros(d-1,1),geodesic);  % Fit a circle at Tp
    if r > pi
        %disp(['Message from getSubSphere.m: (svd) large radius = ' num2str(r) ]);
        r = pi/2;
        [U, dd] = svd(TpData);
        newCenterTp = U(:,end)*pi/2;
    end
    newCenter = ExpNPd(newCenterTp); % Bring back to the sphere by Exp map
    center = rot\newCenter;  % (inv(rot)*newCenter) rotate back the newCenter
    Gnext = objfn(center,r,data);
    err = abs(Gnow - Gnext) ;
    Gnow = Gnext;
    c0 = center;
    cnt = cnt+1;
    if cnt > 30
       %disp(['Message from getSubSphere.m: (svd) iteration reached 30th with error ' num2str(err)]);
       break;
    end
end
i1save.Gnow = Gnow;
i1save.center= center;
i1save.r = r;

% % last eigenvector of covariance matrix
U = pca(data');
% [U, ~, ~] = svd(cov(data'));
initialCenter = U(:,end);

c0 = initialCenter;
TOL = 1e-10;
cnt =0;
err = 1;
[d, n]= size(data);
Gnow = 1e+10;

while err > TOL
    c0 = c0/norm(c0); % normalize the new candidate
    rot = rotMat(c0);  % rotation matrix : c0 -> North Pole
    TpData = LogNPd(rot*data); % Tangent projection by Log map 
    [newCenterTp, r]= LMFsphereFit(TpData,zeros(d-1,1),geodesic);  % Fit a circle at Tp
    if r > pi
        %disp(['Message from getSubSphere.m: (princomp) large radius = ' num2str(r) ]);
        r = pi/2;
        [U, dd] = svd(TpData);
        newCenterTp = U(:,end)*pi/2;
    end
    newCenter = ExpNPd(newCenterTp); % Bring back to the sphere by Exp map
    center = rot\newCenter;  % (inv(rot)*newCenter) rotate back the newCenter
    Gnext = objfn(center,r,data);
    err = abs(Gnow - Gnext) ;
    Gnow = Gnext;
    c0 = center;
    cnt = cnt+1;
    if cnt > 30
       %disp(['Message from getSubSphere.m: (princomp) iteration reached 30th with error ' num2str(err)]);
       break;
    end
end

if i1save.Gnow == min(Gnow,i1save.Gnow)
    center = i1save.center;
    r = i1save.r;
% elseif i2save.Gnow == min([ Gnow i1save.Gnow i2save.Gnow]);
%     disp(['circfit wins psc' num2str([ Gnow i1save.Gnow i2save.Gnow])]);
%     center = i2save.center;
%     r = i2save.r;
% else
%     disp(['princomp wins psc' num2str([ Gnow i1save.Gnow i2save.Gnow])]);
end

if r > pi/2
    center = -center; 
    r = pi - r; 
end % adjust radius
        
% the objective function that we want to minimize: sum of squared distances
% from the data to the subsphere
function g =  objfn(center,r,data) 
  g =  mean((acos(center'*data)-r).^2);

% 
% % initial circle fit function (obsolete) 
% function c = initilalc(x) 
%  % x :  d x n data matrix
% xc = bsxfun(@minus,x,mean(x,2));  % centered x
% A = xc*xc';
% x2 = x.^2;
% x2c = bsxfun(@minus,x2,mean(x2,2));
% if rank(A) < size(x,1)
%     c = pinv(A)*sum(xc*(x2c)',2);
% else
% c = A\sum(xc*(x2c)',2);
% end
% % c = c / d;
% c = c / norm(c);

