function out = Phase_Boxplot(out_warp, alpha, k_p, figs)
% PHASE_BOXPLOT Functional phase Boxplot
% -------------------------------------------------------------------------
% This function constructs the phase boxplot
%
% Usage:  out = Phase_Boxplot(out_warp, k_p, alpha, figs)
%
% Input:
% warp_median: structure from time_warping_median of aligned data using the median
% alpha quantile: value (e.g.,=.05, i.e., 95\%)
% kp: scalar for outlier cutoff (e.g.,=1)
% figs: shows plots of functions (e.g., = true)
% 
% Output:
% Returns a structure containing
% median_x: median warping function
% Q1: First quartile
% Q3: Second quartile
% Q1a: First quantile based on alpha
% Q3a: Second quantile based on alpha
% minn: minimum extreme function
% maxx: maximum extreme function
% outlier_index: indexes of outlier functions

gam = out_warp.gam;

[M, N] = size(gam);
t = linspace(0,1,M);
lambda = 0.5;

% compute phase median
[median_x, psi_median, psi] = SqrtMedian(gam);

% compute phase distances
binsize = mean(diff(t));
dx = zeros(1,N);
v = zeros(M,N);
for i = 1:N
    psi(:,i) = sqrt(gradient(gam(:,i),binsize));
    v(:,i) = inv_exp_map(psi_median,psi(:,i));
    dx(i) = sqrt(trapz(t,v(:,i).^2));
end
[~, dx_ordering] = sort(dx);
CR_50 = dx_ordering(1:ceil(N/2));   % 50% Central Region
m = max(dx(CR_50));                 % Maximal phase distance within 50% Central Region

% identify phase quartiles
angle = zeros(length(CR_50), length(CR_50));
energy = zeros(length(CR_50), length(CR_50));
for i = 1:(length(CR_50)-1)
    for j = (i+1):length(CR_50)
        q1 = v(:,CR_50(i));
        q3 = v(:,CR_50(j));
        q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
        q3=q3/sqrt(trapz(t,q3.^2));
        angle(i,j)=trapz(t,q1.*q3);
        energy(i,j) = (1-lambda) * (dx(CR_50(i))/m + dx(CR_50(j))/m) - lambda * (angle(i,j) + 1);
    end
end
[~, maxloc] = max(energy(:));
[maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);

Q1_index = CR_50(maxloc_row);
Q3_index = CR_50(maxloc_col);
Q1 = gam(:, Q1_index);
Q3 = gam(:, Q3_index);
Q1_psi = sqrt(gradient(Q1,1/(M-1)))';
Q3_psi = sqrt(gradient(Q3,1/(M-1)))';

% identify phase quantiles
[~, dx_ordering] = sort(dx);
CR_alpha = dx_ordering(1:round(N*(1-alpha)));   % 50% Central Region
m = max(dx(CR_alpha));                 % Maximal phase distance within 50% Central Regionles

angle = zeros(length(CR_alpha), length(CR_alpha));
energy = zeros(length(CR_alpha), length(CR_alpha));
for i = 1:(length(CR_alpha)-1)
    for j = (i+1):length(CR_alpha)
        q1 = v(:,CR_alpha(i));
        q3 = v(:,CR_alpha(j));
        q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
        q3=q3/sqrt(trapz(t,q3.^2));
        angle(i,j)=trapz(t,q1.*q3);
        energy(i,j) = (1-lambda) * (dx(CR_alpha(i))/m + dx(CR_alpha(j))/m) - lambda * (angle(i,j) + 1);
    end
end
[~, maxloc] = max(energy(:));
[maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);

Q1a_index = CR_alpha(maxloc_row);
Q3a_index = CR_alpha(maxloc_col);
Q1a = gam(:, Q1a_index);
Q3a = gam(:, Q3a_index);
Q1a_psi = sqrt(gradient(Q1a,1/(M-1)))';
Q3a_psi = sqrt(gradient(Q3a,1/(M-1)))';

% check quartile and quatnile going same direction
tst = trapz(t, v(:,Q1a_index).*v(:,Q1_index));
if (tst < 0)
    Q1a = gam(:,Q3a_index);
    Q3a = gam(:,Q1a_index);
end

% compute phase whiskers
IQR = dx(Q1_index) + dx(Q3_index);
v3 = v(:, Q3_index);
v1 = v(:, Q1_index);
upper_v = v3 + k_p *IQR * v3 / sqrt(trapz(t,v3.^2));
lower_v = v1 + k_p *IQR * v1 / sqrt(trapz(t,v1.^2));
lower_psi = exp_map(psi_median, lower_v);
upper_psi = exp_map(psi_median, upper_v);
lower = cumtrapz(t,lower_psi.^2)';
upper = cumtrapz(t,upper_psi.^2)';

upper_dis = sqrt(trapz(t,(upper_v).^2));
lower_dis = sqrt(trapz(t,(lower_v).^2));
whisker_dis = max(lower_v,upper_v);

% identify phase outliers
outlier_index = [];
for i = 1:N
    if dx(dx_ordering(N+1-i)) > whisker_dis
        outlier_index = [outlier_index; dx_ordering(N+1-i)];
    else
        break
    end
end

% identify phase extremes
distance_to_upper=inf(1,N);
distance_to_lower=inf(1,N);
out_50_CR = setdiff(setdiff([1:N], CR_50), outlier_index);
for i = 1:length(out_50_CR)
    j = out_50_CR(i);
    distance_to_upper(j) = sqrt(trapz(t,(upper_v - v(:,j)).^2));
    distance_to_lower(j) = sqrt(trapz(t,(lower_v - v(:,j)).^2));
end
[~, max_index] = min(distance_to_upper);
[~, min_index] = min(distance_to_lower);
minn = gam(:,min_index);
maxx = gam(:,max_index);
min_psi = psi(:,min_index);
max_psi = psi(:,max_index);

s = linspace(0,1,100);
t = t(:);
median_x = median_x(:);
Fs2 = zeros(length(t), 595);
Fs2(:,1) = (1-s(1)) * (minn-t) + s(1) * (Q1-t);
for j=2:100
    Fs2(:,j) = (1-s(j)) * (minn-t) + s(j) * (Q1a-t);
    Fs2(:,99+j) = (1-s(j)) * (Q1a-t) + s(j) * (Q1-t);
    Fs2(:,198+j) = (1-s(j)) * (Q1-t) + s(j) * (median_x-t);
    Fs2(:,297+j) = (1-s(j)) * (median_x-t) + s(j) * (Q3-t);
    Fs2(:,396+j) = (1-s(j)) * (Q3-t) + s(j) * (Q3a-t);
    Fs2(:,495+j) = (1-s(j)) * (Q3a-t) + s(j) * (maxx-t);
end
Q1_psi = Q1_psi(:);
Q1a_psi = Q1a_psi(:);
Q3_psi = Q3_psi(:);
Q3a_psi = Q3a_psi(:);
d1=acos(trapz(t,psi_median.*Q1_psi));
d1a=acos(trapz(t,Q1_psi.*Q1a_psi));
dl=acos(trapz(t,Q1a_psi.*min_psi));
d3=acos(trapz(t,psi_median.*Q3_psi));
d3a=acos(trapz(t,Q3_psi.*Q3a_psi));
du=acos(trapz(t,Q3a_psi.*max_psi));
part1=linspace(-d1-d1a-dl,-d1-d1a,100);
part2=linspace(-d1-d1a,-d1,100);
part3=linspace(-d1,0,100);
part4=linspace(0,d3,100);
part5=linspace(d3,d3+d3a,100);
part6=linspace(d3+d3a,d3+d3a+du,100);
allparts=[part1,part2(2:100),part3(2:100),part4(2:100),part5(2:100),part6(2:100)];
[U,V]=meshgrid(linspace(0,1,M),allparts);
U=U.';
V=V.';

plt.U=U;
plt.V=V;
plt.allparts = allparts;
plt.Fs2 = Fs2;
plt.d1 = d1;
plt.d1a = d1a;
plt.dl = dl;
plt.d3 = d3;
plt.d3a = d3a;
plt.du = du;
plt.Q1_psi = Q1a_psi;
plt.Q3_psi = Q3a_psi;

if (figs)

    figure(410); clf;
    plot(t, median_x, 'black','linewidth', 2);
    hold on;
    plot(t, Q1, 'blue','linewidth', 2);
    plot(t, Q3, 'blue', 'linewidth', 2);
    plot(t, Q1a, 'green', 'linewidth', 2);
    plot(t, Q3a, 'green', 'linewidth', 2);
    plot(t,maxx,'red','linewidth',2);
    plot(t,minn,'red','linewidth',2);
    axis square;
    axis([0,1,0,1]);

    figure(416); clf;
    surf(U,V,Fs2);
    hold on;
    shading flat;
    plot3(t,zeros(1,M),median_x - t,'k','LineWidth',3)
    plot3(t,repmat(-d1,M,1),Q1 - t,'b','LineWidth',3)
    plot3(t,repmat(-d1-d1a,M,1),Q1a - t,'g','LineWidth',3)
    plot3(t,repmat(-d1-d1a-dl,M,1),minn -t,'r','LineWidth',3)
    plot3(t,repmat(d3,M,1),Q3 - t,'b','LineWidth',3)
    plot3(t,repmat(d3+d3a,M,1),Q3a - t,'g','LineWidth',3)
    plot3(t,repmat(d3+d3a+du,M,1),maxx - t,'r','LineWidth',3)
    axis square;
end

out.time = time;
out.psi_median = psi_median;
out.Q1 = Q1;
out.Q3 = Q3;
out.Q1a = Q1a;
out.minn = minn;
out.maxx = maxx;
out.outlier_index = outlier_index;
out.plt = plt;
out.median_x = median_x;
