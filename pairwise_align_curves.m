function [dist,X2n,q2n,X1,q1]=pairwise_align_curves(X1,X2,option)
% PAIRWISE_ALIGN_CURVES registers two curves
% -------------------------------------------------------------------------
% Convert to SRSF
%
% Usage: q = f_to_srvf(f,time)
%
% This function converts functions to srsf
%
% Input:
% X1: matrix (nxT) of a n dimensional curve of T points
% X2: matrix (nxT) of a n dimensional curve of T points
%
% default options
% option.plot_geod = true;  % plot geodesic
% option.plot_reg = true;   % plot registration
% option.plot_reparam = true;  % plot reparametization
% option.print = true;  % print output
% option.N = 200;  % resample curve to N points
% option.stp = 6;  % number of steps in geodesic
% option.closed = false;  % closed curve
%
% Output:
% q: matrix of SRSFs(X1,X2)
% dist: distance between two curves
% X2n: optimally reigstered curve X2
% q2n: optimally registered SRVF
% X1: normalized curve X1
% q1: SRVF of of curve 1

% Load some parameters, no need to change this
n = size(X1,1);

if nargin < 3
    option.plot_geod = true;
    option.plot_reg = true;
    option.plot_reparam = true;
    option.print = true;
    option.N = 200;
    option.stp = 6;
    option.closed = false;
end

% Resample the curves to have N points
X1 = ReSampleCurve(X1,option.N);
X2 = ReSampleCurve(X2,option.N);

%Center curves, not really needed but good for display purposes
a = -calculateCentroid(X1);
X1 = X1 + repmat(a,1,option.N);
a = -calculateCentroid(X2);
X2 = X2 + repmat(a,1,option.N);

% Form the q function for representing curves and find best rotation
q1 = curve_to_q(X1);
q2 = curve_to_q(X2);
[~,R,gamI] = Find_Rotation_and_Seed_unique(q1,q2,true,option.closed);
X2 = R*X2;
X2n = warp_curve_gamma(X2,gamI);
q2n = curve_to_q(X2n);

if option.plot_reparam
    figure(100); clf;
    plot(linspace(0,1,option.N),gamI,'LineWidth',2)
    axis square
    grid on
    title('Warping function')
end

% Find optimal rotation
[q2n,R] = Find_Best_Rotation(q1,q2n);
X2n = R*X2n;

% Forming geodesic between the registered curves
q1dotq2=InnerProd_Q(q1,q2n);

% Compute shooting vector
if q1dotq2>1
    q1dotq2=1;
end

dist = acos(q1dotq2);

if option.print
    fprintf('The distance between the two curves is %0.3f\n',dist)
end
if(option.plot_geod)
    PsiQ = geodesic_sphere_Full(q1,q2,option.stp,option.closed);
    p2n = q_to_curve(q2n);
    Path_Plot(PsiQ,p2n,10,[73,6])
    Geod = zeros(n,option.N,option.stp+1);
    for j=1:option.stp+1
        Geod(:,:,j)=q_to_curve(PsiQ(:,:,j));
    end
end
% Displaying the correspondence
if(option.plot_reg)
    X2n=Geod(:,:,end);
    X1=Geod(:,:,1);
    if (size(X1,1)==2)
        figure(3); clf;
        z = plot(X1(1,:), X1(2,:),'r');
        set(z,'LineWidth',2);
        axis off;
        hold on;
        z = plot(X2n(1,:), 0.15+X2n(2,:),'b-+');
        set(z,'LineWidth',3);
        N = size(X1,2);
        for j=1:N/15
            i = j*15;
            plot([X1(1,i) X2n(1,i)],[X1(2,i) 0.15+X2n(2,i)], 'k');
        end
    else
        figure(3); clf;
        z = plot3(X1(1,:), X1(2,:), X1(3,:),'r');
        set(z,'LineWidth',2);
        axis off;
        hold on;
        z = plot3(X2n(1,:), 0.15+X2n(2,:), X2n(3,:),'b-+');
        set(z,'LineWidth',3);
        N = size(X1,2);
        for j=1:N/15
            i = j*15;
            plot3([X1(1,i) X2n(1,i)],[X1(2,i) 0.15+X2n(2,i)], [X1(3,i) X2n(3,i)], 'k');
        end
    end
end