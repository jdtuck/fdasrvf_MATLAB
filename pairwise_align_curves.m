function [dist,X2n,q2n,X1,q1]=pairwise_align_curves(X1,X2,option)
% PAIRWISE_ALIGN_CURVES registers to curves
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
% option.N = 200;  % resample curve to N points
%
% Output:
% q: matrix of SRSFs(X1,X2)
% dist: distance between two curves
% X2n: optimally reigstered curve X2
% q2n: optimally registered SRVF
% X1: normalized curve X1
% q1: SRVF of of curve 1

% Load some parameters, no need to change this
lam = 0;
n = size(X1,1);

if nargin < 3
    option.plot_geod = true;
    option.plot_reg = true;
    option.plot_reparam = true;
    option.N = 200;
end

% Resample the curves to have N points
X1 = ReSampleCurve(X1,option.N);
X2 = ReSampleCurve(X2,option.N);

%Center curves, not really needed but good for display purposes
X1 = X1 - repmat(mean(X1,2),1,size(X1,2));
X2 = X2 - repmat(mean(X2,2),1,size(X2,2));

% Form the q function for representing curves and find best rotation
q1 = curve_to_q(X1);
q2 = curve_to_q(X2);
A = q1*q2';
[U,~,V] = svd(A);
if det(A)> 0
    Ot = U*V';
else
    if (size(X1,1)==2)
        Ot = U*([V(:,1) -V(:,2)])';
    else
        Ot = U*([V(:,1) V(:,2) -V(:,3)])';
    end
end
X2 = Ot*X2;
q2 = Ot*q2;

% Applying optimal re-parameterization to the second curve
[gam] = DynamicProgrammingQ(q1/sqrt(InnerProd_Q(q1,q1)),q2/sqrt(InnerProd_Q(q2,q2)),lam,0);
gamI = invertGamma(gam);
gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
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
A = q1*q2n';
[U,~,V] = svd(A);
if det(A)> 0
    Ot = U*V';
else
    if (size(X1,1)==2)
        Ot = U*([V(:,1) -V(:,2)])';
    else
        Ot = U*([V(:,1) V(:,2) -V(:,3)])';
    end
end
X2n = Ot*X2n;
q2n = Ot*q2n;

% Forming geodesic between the registered curves
N = size(X1,2);
dist = acos(sum(sum(q1.*q2n))/N);
fprintf('The distance between the two curves is %0.3f\n',dist)
no = 7;
PsiQ = zeros(n,option.N,no);
PsiX = zeros(n,option.N,no);
if(option.plot_geod)
    if (size(X1,1)==2)
        for t=1:no
            s = dist*(t-1)/6;
            PsiQ(:,:,t) = (sin(dist - s)*q1 + sin(s)*q2n)/sin(dist);
            PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
        end
        figure(4); clf; axis equal; hold on;
        for t=1:7
            z = plot(0.5*t + PsiX(1,:,t), PsiX(2,:,t),'r-');
            set(z,'LineWidth',2,'color',[(t-1)/6 (t-1)/12 0]);
        end
        axis off;
    else
        for t=1:7
            s = dist*(t-1)/6;
            PsiQ(:,:,t) = (sin(dist - s)*q1 + sin(s)*q2n)/sin(dist);
            PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
        end
        figure(4); clf; axis equal; hold on;
        for t=1:7
            z = plot3(0.5*t + PsiX(1,:,t), PsiX(2,:,t), PsiX(3,:,t),'r-');
            set(z,'LineWidth',2,'color',[(t-1)/6 (t-1)/12 0]);
        end
        axis off;
    end
end

% Displaying the correspondence
if(option.plot_reg)
    X2n=PsiX(:,:,end);
    X1=PsiX(:,:,1);
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