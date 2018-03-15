function out = Amplitude_Boxplot(out_warp, alpha, k_a, figs)
% AMPLITUDE_BOXPLOT Functional Amplitude Boxplot
% -------------------------------------------------------------------------
%
% This function constructs the amplitude boxplot
%
% Usage:  out = Amplitude_Boxplot(out_warp, alpha, k_a, figs)
%
% Input:
% out_warp: struct from time_warping_median of aligned data using the median
% alpha: quantile value (e.g.,=.05, i.e., 95\%)
% ka: scalar for outlier cutoff (e.g.,=1)
% figs: shows plots of functions (e.g., = true)
%
% Output: structure containing
% median_y: median function
% Q1: First quartile
% Q3: Second quartile
% Q1a: First quantile based on alpha
% Q3a: Second quantile based on alpha
% minn: minimum extreme function
% maxx: maximum extreme function
% outlier_index: indexes of outlier functions
% fmedian: median function

if out_warp.rsamps
    f_tilde = out_warp.fs;
    f_median = out_warp.fmedian;
    q_tilde = out_warp.qs;
    q_median = out_warp.mqn;
    t = out_warp.time;
else
    f_tilde = out_warp.fn;
    f_median = out_warp.fmedian;
    q_tilde = out_warp.qn;
    q_median = out_warp.mqn;
    t = out_warp.time; 
end

[M, N] = size(f_tilde);
lambda = 0.5;

% amplitude median
median_y = f_median;

% compute amplitude distances
dy = zeros(1,N);
for i = 1:N
    dy(i) = sqrt(trapz(t,(q_median-q_tilde(:,i)).^2));
end
[~, dy_ordering] = sort(dy);
CR_50 = dy_ordering(1:ceil(N/2));       % 50% Central Region
m = max(dy(CR_50));                     % Maximal amplitude distance within 50% Central Region

% identify amplitude quartiles
angle = zeros(length(CR_50), length(CR_50));
energy = zeros(length(CR_50), length(CR_50));
for i = 1:(length(CR_50)-1)
    for j = (i+1):length(CR_50)
        q1 = q_tilde(:,CR_50(i)) - q_median;
        q3 = q_tilde(:,CR_50(j)) - q_median;
        q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
        q3=q3/sqrt(trapz(t,q3.^2));
        angle(i,j)=trapz(t,q1.*q3);
        energy(i,j) = (1-lambda) * (dy(CR_50(i))/m + dy(CR_50(j))/m) - lambda * (angle(i,j) + 1);
    end
end
[~, maxloc] = max(energy(:));
[maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);

Q1_index = CR_50(maxloc_row);
Q3_index = CR_50(maxloc_col);
Q1_q = q_tilde(:,Q1_index);
Q3_q = q_tilde(:,Q3_index);
Q1 = f_tilde(:,Q1_index);
Q3 = f_tilde(:,Q3_index);

% identify amplitude quantiles
[~, dy_ordering] = sort(dy);
CR_alpha = dy_ordering(1:round(N*(1-alpha)));       % (1-alpha)% Central Region
m = max(dy(CR_alpha));                     % Maximal amplitude distance within (1-alpha)% Central Region
angle = zeros(length(CR_alpha), length(CR_alpha));
energy = zeros(length(CR_alpha), length(CR_alpha));
for i = 1:(length(CR_alpha)-1)
    for j = (i+1):length(CR_alpha)
        q1 = q_tilde(:,CR_alpha(i)) - q_median;
        q3 = q_tilde(:,CR_alpha(j)) - q_median;
        q1=q1/sqrt(trapz(t,q1.^2));     % normalize to unit 1
        q3=q3/sqrt(trapz(t,q3.^2));
        angle(i,j)=trapz(t,q1.*q3);
        energy(i,j) = (1-lambda) * (dy(CR_alpha(i))/m + dy(CR_alpha(j))/m) - lambda * (angle(i,j) + 1);
    end
end
[~, maxloc] = max(energy(:));
[maxloc_row, maxloc_col] = ind2sub(size(energy), maxloc);

Q1a_index = CR_alpha(maxloc_row);
Q3a_index = CR_alpha(maxloc_col);
Q1a_q = q_tilde(:,Q1a_index);
Q3a_q = q_tilde(:,Q3a_index);
Q1a = f_tilde(:,Q1a_index);
Q3a = f_tilde(:,Q3a_index);

% compute amplitude whiskers
IQR = dy(Q1_index)+dy(Q3_index);
v1 = Q1_q - q_median;
v3 = Q3_q - q_median;
upper_q = Q3_q + k_a * IQR * v3 / sqrt(trapz(t,v3.^2));
lower_q = Q1_q + k_a * IQR * v1 / sqrt(trapz(t,v1.^2));

upper_dis = sqrt(trapz(t,(upper_q - q_median).^2));
lower_dis = sqrt(trapz(t,(lower_q - q_median).^2));
whisker_dis = max([lower_dis upper_dis]);

% identify amplitude outliers
outlier_index = [];
for i = 1:N
    if dy(dy_ordering(N+1-i)) > whisker_dis
        outlier_index = [outlier_index; dy_ordering(N+1-i)];
    else
        break
    end
end

% identify amplitude extremes
distance_to_upper=inf(1,N);
distance_to_lower=inf(1,N);
out_50_CR = setdiff(setdiff((1:N), CR_50), outlier_index);
for i = 1:length(out_50_CR)
    j = out_50_CR(i);
    distance_to_upper(j) = sqrt(trapz(t,(upper_q - q_tilde(:,j)).^2));
    distance_to_lower(j) = sqrt(trapz(t,(lower_q - q_tilde(:,j)).^2));
end
[~, max_index] = min(distance_to_upper);
[~, min_index] = min(distance_to_lower);
min_q = q_tilde(:,min_index);
max_q = q_tilde(:,max_index);
minn = f_tilde(:,min_index);
maxx = f_tilde(:,max_index);

s = linspace(0,1,100);
Fs2 = zeros(length(t), 595);
Fs2(:,1) = (1-s(1)) * minn + s(1) * Q1;    % Final surface plot
for j=2:100
    Fs2(:,j) = (1-s(j)) * minn + s(j) * Q1a;
    Fs2(:,99+j) = (1-s(j)) * Q1a + s(j) * Q1;
    Fs2(:,198+j) = (1-s(j)) * Q1 + s(j) * f_median;
    Fs2(:,297+j) = (1-s(j)) * f_median + s(j) * Q3;
    Fs2(:,396+j) = (1-s(j)) * Q3 + s(j) * Q3a;
    Fs2(:,495+j) = (1-s(j)) * Q3a + s(j) * maxx;
end
d1=sqrt(trapz(t,(q_median-Q1_q).^2));
d1a=sqrt(trapz(t,(Q1_q-Q1a_q).^2));
dl=sqrt(trapz(t,(Q1a_q-min_q).^2));
d3=sqrt(trapz(t,(q_median-Q3_q).^2));
d3a=sqrt(trapz(t,(Q3_q-Q3a_q).^2));
du=sqrt(trapz(t,(Q3a_q-max_q).^2));
part1=linspace(-d1-d1a-dl,-d1-d1a,100);
part2=linspace(-d1-d1a,-d1,100);
part3=linspace(-d1,0,100);
part4=linspace(0,d3,100);
part5=linspace(d3,d3+d3a,100);
part6=linspace(d3+d3a,d3+d3a+du,100);
allparts=[part1,part2(2:100),part3(2:100),part4(2:100),part5(2:100),part6(2:100)];
[U,V]=meshgrid(t,allparts);
U=U';
V=V';

plt.U=U;
plt.V=V;
plt.Fs2 = Fs2;
plt.allparts = allparts;
plt.d1 = d1;
plt.d1a = d1a;
plt.dl = dl;
plt.d3 = d3;
plt.d3a = d3a;
plt.du = du;
plt.Q1q = Q1a_q;
plt.Q3q = Q3a_q;

if (figs)
    figure(310); clf;
    plot(t, f_median, 'black','linewidth', 2);
    hold on;
    plot(t, Q1, 'blue','linewidth', 2);
    plot(t, Q3, 'blue', 'linewidth', 2);
    plot(t, Q1a, 'green', 'linewidth', 2);
    plot(t, Q3a, 'green', 'linewidth', 2);
    plot(t,minn,'red', 'linewidth',2);
    plot(t,maxx,'red', 'linewidth',2);
    xlim([t(1) t(end)]);
    ylim auto;
    
    
    figure(311); clf;
    surf(U,V,Fs2);
    hold on;
    shading flat;
    plot3(t,zeros(1,M),f_median,'k','LineWidth',3)
    plot3(t,repmat(-d1,M,1),Q1,'b','LineWidth',3)
    plot3(t,repmat(-d1-d1a,M,1),Q1a,'g','LineWidth',3)
    plot3(t,repmat(-d1-d1a-dl,M,1),minn,'r','LineWidth',3)
    plot3(t,repmat(d3,M,1),Q3,'b','LineWidth',3)
    plot3(t,repmat(d3+d3a,M,1),Q3a,'g','LineWidth',3)
    plot3(t,repmat(d3+d3a+du,M,1),maxx,'r','LineWidth',3)
end

out.f_median = median_y;
out.q_median = q_median;
out.Q1 = Q1;
out.Q3 = Q3;
out.Q1a = Q1a;
out.Q3a = Q3a;
out.minn = minn;
out.maxx = maxx;
out.outlier_index = outlier_index;
out.plt = plt;
