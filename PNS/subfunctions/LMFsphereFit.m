function [center, r] = LMFsphereFit(A, initialCenter, geodesic)
% LMFSPHEREFIT : The least square estimates of the sphere to the data.
%              The non-linear least squares problem is solved with the
%              Levenberg-Marquardt algorithm of the Optimization Toolbox.
% [center]= LMFsphereFit(A) with d x n data matrix A (any d = 2, 3, ...).
% [center, r]= LMFsphereFit(A) with d x n data matrix gives the center and
%                              the radius.
% [center, r]= LMFsphereFit(A,initialCenter,1) forces the sphere radius to
%                              be pi/2.
%
% % example;
% n =50;
% theta = linspace(0,pi*1.5,n);
% data = 5*[cos(theta); sin(theta)] + randn(2,n);
% [x,r]=LMFsphereFit(data);
% figure(1);clf
% scatter(data(1,:),data(2,:),'.b');hold on;
% plot(x(1),x(2),'or');
% Estcirc = r*[cos(theta); sin(theta)]+repmat(x,1,size(data,2));
% plot(Estcirc(1,:),Estcirc(2,:))
% axis equal
%
%   See also lsqnonlin.

% Last updated Aug 10, 2009
% Sungkyu Jung

if nargin < 2 || isempty(initialCenter)
    initialCenter = mean(A,2);
end
if nargin < 3
    geodesic = 0; % Do not force the radius to be pi/2
end

opts = optimoptions("lsqnonlin", Algorithm="levenberg-marquardt", ...
    Display="off", MaxIterations=50, StepTolerance=1e-9);
center = lsqnonlin(@(c) sphereRes(c,A,geodesic), initialCenter, [], [], opts);

r = sphereRadius(center, A, geodesic);
end

function rr = sphereRes(center, A, geodesic)
% residuals of the sphere fit for the given center
di = vecnorm(A - center, 2, 1);
rr = (di - sphereRadius(center, A, geodesic)).';
end

function r = sphereRadius(center, A, geodesic)
if geodesic == 0
    r = mean(vecnorm(A - center, 2, 1));
else
    r = pi/2;
end
end
