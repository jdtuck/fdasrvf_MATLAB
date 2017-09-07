function [median_x, Q1, Q3, minn, maxx, outlier_index] = Phase_Boxplot(gam, k_p)
% k_p: outlier cutoff constant for phase component

[N, M] = size(gam);
t = linspace(0,1,M);
lambda = 0.5;
figs = 0;

if figs==1
    figure(400); clf;
        plot(linspace(0,1,M), gam', 'linewidth', 2);
        axis square;
        axis([0,1,0,1]);
        set(gca,'FontSize',19);
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
end

% compute phase median
[median_x, psi_median, psi] = phase_median(gam);

% compute phase distances
binsize = mean(diff(t));
for i = 1:N
    psi(:,i) = sqrt(gradient(gam(i,:),binsize));
    v(:,i) = inv_exp_map(psi_median,psi(:,i));
    dx(i) = sqrt(trapz(t,v(:,i).^2));
end
[~, dx_ordering] = sort(dx);
CR_50 = dx_ordering(1:ceil(N/2));   % 50% Central Region
m = max(dx(CR_50));                 % Maximal phase distance within 50% Central Region

% identify phase quartiles
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
[maxval, maxloc] = max(energy(:));
[maxloc_row maxloc_col] = ind2sub(size(energy), maxloc);

Q1_index = CR_50(maxloc_row);
Q3_index = CR_50(maxloc_col);
Q1 = gam(Q1_index,:);
Q3 = gam(Q3_index,:);
Q1_psi = sqrt(gradient(Q1,1/(M-1)))';
Q3_psi = sqrt(gradient(Q3,1/(M-1)))';

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

if figs == 1
    figure(401); clf;
        plot(t, median_x, 'black','linewidth', 2);
        axis square;
        axis([0,1,0,1]);
        hold on;
        set(gca,'FontSize',19)
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
        plot(t, Q1, 'blue','linewidth', 2); 
        plot(t, Q3,'green', 'linewidth', 2);
        plot(t,upper,'magenta','linewidth',2)
        plot(t,lower,'red','linewidth',2)

    % surface plot
    s = linspace(0,1,100);
    Fs1(:,1) = (1-s(1)) * (lower-t) + s(1) * (Q1-t); 
    for j=2:100
        Fs1(:,j) = (1-s(j)) * (lower-t) + s(j) * (Q1-t); 
        Fs1(:,99+j) = (1-s(j)) * (Q1-t) + s(j) * (median_x-t); 
        Fs1(:,198+j) = (1-s(j)) * (median_x-t) + s(j) * (Q3-t); 
        Fs1(:,297+j) = (1-s(j)) * (Q3-t) + s(j) * (upper-t); 
    end
    dl=acos(trapz(t,lower_psi.*Q1_psi));
    d1=acos(trapz(t,Q1_psi.*psi_median));
    d3=acos(trapz(t,psi_median.*Q3_psi));
    du=acos(trapz(t,Q3_psi.*upper_psi));
    part1=linspace(-d1-dl,-d1,100);
    part2=linspace(-d1,0,100);
    part3=linspace(0,d3,100);
    part4=linspace(d3,d3+du,100);
    allparts=[part1,part2(2:100),part3(2:100),part4(2:100)];
    [U,V]=meshgrid(t,allparts);
    
    figure(402); clf;
        surf(U',V',Fs1);
        hold on;
        shading flat;
        light;
        plot3(t,zeros(1,M),median_x - t,'k','LineWidth',3);
        plot3(t,repmat(-d1,M,1),Q1 - t,'b','LineWidth',3);
        plot3(t,repmat(-d1-dl,M,1),lower -t,'r','LineWidth',3);
        plot3(t,repmat(d3,M,1),Q3 - t,'g','LineWidth',3);
        plot3(t,repmat(d3+du,M,1),upper - t,'magenta','LineWidth',3);
        set(gca,'FontSize',19)
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
end

upper_dis = sqrt(trapz(t,(upper_v).^2));
lower_dis = sqrt(trapz(t,(lower_v).^2));
whisker_dis = max(lower_v,upper_v);

% identify phase outliers
outlier_index = [];
for i = 1:N
    if dx(dx_ordering(N+1-i)) > whisker_dis
        outlier_index = [outlier_index; dx_ordering(N+1-i)];
    else break
    end
end

% identify phase extremes
distance_to_upper=repmat(inf,1,N);
distance_to_lower=repmat(inf,1,N);
out_50_CR = setdiff(setdiff([1:N], CR_50), outlier_index);
for i = 1:length(out_50_CR)
    j = out_50_CR(i);
    distance_to_upper(j) = sqrt(trapz(t,(upper_v - v(:,j)).^2)); 
    distance_to_lower(j) = sqrt(trapz(t,(lower_v - v(:,j)).^2)); 
end;
[~, max_index] = min(distance_to_upper);
[~, min_index] = min(distance_to_lower);
minn = gam(min_index,:);
maxx = gam(max_index,:);
min_psi = psi(:,min_index);
max_psi = psi(:,max_index);

figure(410); clf; 
    plot(t, median_x, 'black','linewidth', 2);
    hold on;
    plot(t, Q1, 'blue','linewidth', 2); 
    plot(t, Q3,'green', 'linewidth', 2);
    plot(t,maxx,'magenta','linewidth',2);
    plot(t,minn,'red','linewidth',2);
    axis square;
    axis([0,1,0,1]);
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

s = linspace(0,1,100);
Fs2(:,1) = (1-s(1)) * (minn-t) + s(1) * (Q1-t); 
for j=2:100
    Fs2(:,j) = (1-s(j)) * (minn-t) + s(j) * (Q1-t); 
    Fs2(:,99+j) = (1-s(j)) * (Q1-t) + s(j) * (median_x-t); 
    Fs2(:,198+j) = (1-s(j)) * (median_x-t) + s(j) * (Q3-t); 
    Fs2(:,297+j) = (1-s(j)) * (Q3-t) + s(j) * (maxx-t); 
end
dl=acos(trapz(t,min_psi.*Q1_psi));
d1=acos(trapz(t,Q1_psi.*psi_median));
d3=acos(trapz(t,psi_median.*Q3_psi));
du=acos(trapz(t,Q3_psi.*max_psi));
part1=linspace(-d1-dl,-d1,100);
part2=linspace(-d1,0,100);
part3=linspace(0,d3,100);
part4=linspace(d3,d3+du,100);
allparts=[part1,part2(2:100),part3(2:100),part4(2:100)];
[U,V]=meshgrid(linspace(0,1,M),allparts);

figure(416); clf;
    surf(U',V',Fs2);
    hold on;
    shading flat;
    light;
    plot3(t,zeros(1,M),median_x - t,'k','LineWidth',3)
    plot3(t,repmat(-d1,M,1),Q1 - t,'b','LineWidth',3)
    plot3(t,repmat(-d1-dl,M,1),minn -t,'r','LineWidth',3)
    plot3(t,repmat(d3,M,1),Q3 - t,'g','LineWidth',3)
    plot3(t,repmat(d3+du,M,1),maxx - t,'magenta','LineWidth',3)
    axis square;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

for i = 1:length(outlier_index)   % display phase outliers
    figure(420 + i); clf;
        for j = 1:N
            plot(t,gam(j,:),'green','linewidth',2);
            hold on;
        end
        axis square;
        axis([0,1,0,1]);
        set(gca,'FontSize',19)
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
        plot(t,median_x,'black','linewidth',2);
        plot(t,gam(outlier_index(i),:),'red','linewidth',2);
end