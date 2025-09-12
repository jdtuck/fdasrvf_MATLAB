function [center, r] = LMFsphereFit(A,initialCenter,geodesic)
% LMFSPHEREFIT : The least square estimates of the sphere to the data.
%              The non-linear least square estimates are calculated by
%              the Levenberg-Marquardt method in Fletcher's modification
%              (Fletcher, R., (1971): A Modified Marquardt Subroutine for
%               Nonlinear Least Squares. Rpt. AERE-R 6799, Harwell)
%              and implemented for MATLAB by M. Balda's "LMFnlsq.m"
% [center]= LMFsphereFit(A) with d x n data matrix A (any d = 2, 3, ...).
% [center, r]= LMFsphereFit(A) with d x n data matrix gives the center and 
%                              the radius.
% [center, r]= LMFsphereFit(A,1) forces the sphere radius as 1
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
%   See also LMFcircRes,LMFnlsq.

% Last updated Aug 10, 2009
% Sungkyu Jung

if nargin==0 && nargout==0, help LMFsphereFit, return, end     %   Display help

if nargin==1
    initialCenter =  mean(A,2);
    geodesic = 0; % Do not force the radius to be pi/2
end

if nargin==2
    geodesic = 0; % Do not force the radius to be pi/2
end

global data ; data = A;
global greatCircle; greatCircle = geodesic;

center=LMFnlsq('LMFsphereRes',initialCenter,'Display',0,'MaxIter',50,'XTol',1e-9);
di = sqrt(sum(...
    (data - repmat(center,1,size(data,2))).^2 ...
    ));
clear data;
if geodesic == 0
    r= mean(di);
else 
    r = pi/2;
end
