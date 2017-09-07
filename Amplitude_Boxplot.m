function [median_y, Q1, Q3, minn, maxx, outlier_index] = Amplitude_Boxplot(f_tilde, f_median, q_tilde, q_median, t, k_a)
% k_a: outlier cutoff constant for amplitude component

[M, N] = size(f_tilde);
lambda = 0.5;
figs = 0;

if figs==1
    figure(300); clf;
        plot(t, f_tilde, 'LineWidth',2);
        xlim([t(1) t(end)]);
        ylim auto;
        set(gca,'FontSize',19);
        figax=axis;
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
end

% amplitude median
median_y = f_median;

% compute amplitude distances
for i = 1:N
    dy(i) = sqrt(trapz(t,(q_median-q_tilde(:,i)).^2));
end
[~, dy_ordering] = sort(dy);
CR_50 = dy_ordering(1:ceil(N/2));       % 50% Central Region
m = max(dy(CR_50));                     % Maximal amplitude distance within 50% Central Region

% identify amplitude quartiles
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
[maxval, maxloc] = max(energy(:));
[maxloc_row maxloc_col] = ind2sub(size(energy), maxloc);

Q1_index = CR_50(maxloc_row);
Q3_index = CR_50(maxloc_col);
Q1_q = q_tilde(:,Q1_index);
Q3_q = q_tilde(:,Q3_index);
Q1 = f_tilde(:,Q1_index);
Q3 = f_tilde(:,Q3_index);

% compute amplitude whiskers
IQR = dy(Q1_index)+dy(Q3_index);
v1 = Q1_q - q_median;
v3 = Q3_q - q_median;
upper_q = Q3_q + k_a * IQR * v3 / sqrt(trapz(t,v3.^2));
lower_q = Q1_q + k_a * IQR * v1 / sqrt(trapz(t,v1.^2));
upper = cumtrapz(t,upper_q.*abs(upper_q));
lower = cumtrapz(t,lower_q.*abs(lower_q));

if figs == 1
    figure(301); clf; 
        plot(t, f_median, 'black','linewidth', 2); 
        hold on;
        plot(t, Q1, 'blue','linewidth', 2); 
        plot(t, Q3,'green', 'linewidth', 2);
        plot(t,lower,'red','linewidth',2);
        plot(t,upper,'magenta','linewidth',2)
        xlim([t(1) t(end)]);
        ylim auto;
        set(gca,'FontSize',18);
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]); 
                 
    % surface plot
    s = linspace(0,1,100);
    Fs1(:,1) = (1-s(1)) * lower + s(1) * Q1; 
    for j=2:100
        Fs1(:,j) = (1-s(j)) * lower + s(j) * Q1; 
        Fs1(:,99+j) = (1-s(j)) * Q1 + s(j) * f_median; 
        Fs1(:,198+j) = (1-s(j)) * f_median + s(j) * Q3; 
        Fs1(:,297+j) = (1-s(j)) * Q3 + s(j) * upper; 
    end
    d1=sqrt(trapz(t,(q_median-Q1_q).^2));
    dl=sqrt(trapz(t,(Q1_q-lower_q).^2));
    d3=sqrt(trapz(t,(q_median-Q3_q).^2));
    du=sqrt(trapz(t,(Q3_q-upper_q).^2));
    part1=linspace(-d1-dl,-d1,100);
    part2=linspace(-d1,0,100);
    part3=linspace(0,d3,100);
    part4=linspace(d3,d3+du,100);
    allparts=[part1,part2(2:100),part3(2:100),part4(2:100)];
    [U,V]=meshgrid(t,allparts);

    figure(302); clf;
        surf(U',V',Fs1);
        hold on;
        shading flat;
        light;    
        plot3(linspace(0,1,M),zeros(1,M),f_median,'k','LineWidth',3)
        plot3(linspace(0,1,M),repmat(-d1,M,1),Q1,'b','LineWidth',3)
        plot3(linspace(0,1,M),repmat(-d1-dl,M,1),lower,'r','LineWidth',3)
        plot3(linspace(0,1,M),repmat(d3,M,1),Q3,'g','LineWidth',3)
        plot3(linspace(0,1,M),repmat(d3+du,M,1),upper,'magenta','LineWidth',3)
        set(gca,'FontSize',18);
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
end

upper_dis = sqrt(trapz(t,(upper_q - q_median).^2));
lower_dis = sqrt(trapz(t,(lower_q - q_median).^2));
whisker_dis = max([lower_dis upper_dis]);

% identify amplitude outliers
outlier_index = [];
for i = 1:N
    if dy(dy_ordering(N+1-i)) > whisker_dis
        outlier_index = [outlier_index; dy_ordering(N+1-i)];
    else break
    end
end

% identify amplitude extremes
distance_to_upper=repmat(inf,1,N);
distance_to_lower=repmat(inf,1,N);
out_50_CR = setdiff(setdiff([1:N], CR_50), outlier_index);
for i = 1:length(out_50_CR)
    j = out_50_CR(i);
    distance_to_upper(j) = sqrt(trapz(t,(upper_q - q_tilde(:,j)).^2)); 
    distance_to_lower(j) = sqrt(trapz(t,(lower_q - q_tilde(:,j)).^2)); 
end;
[~, max_index] = min(distance_to_upper);
[~, min_index] = min(distance_to_lower);
min_q = q_tilde(:,min_index);
max_q = q_tilde(:,max_index);
minn = f_tilde(:,min_index);
maxx = f_tilde(:,max_index);

figure(310); clf; 
    plot(t, f_median, 'black','linewidth', 2) 
    hold on;
    plot(t, Q1, 'blue','linewidth', 2) 
    plot(t, Q3,'green', 'linewidth', 2);
    plot(t,minn,'red','linewidth',2);
    plot(t,maxx,'magenta','linewidth',2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',18);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

figure(311); clf; 
    plot(t, minn,'red','linewidth',2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
figure(312); clf; 
    plot(t, Q1, 'blue','linewidth', 2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
figure(313); clf; 
    plot(t, f_median, 'black','linewidth', 2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
figure(314); clf; 
    plot(t, Q3,'green', 'linewidth', 2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
figure(315); clf; 
    plot(t, maxx,'magenta','linewidth',2);
    xlim([t(1) t(end)]);
    ylim auto;
    set(gca,'FontSize',19)
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

s = linspace(0,1,100);
Fs2(:,1) = (1-s(1)) * minn + s(1) * Q1;    % Final surface plot
for j=2:100
    Fs2(:,j) = (1-s(j)) * minn + s(j) * Q1; 
    Fs2(:,99+j) = (1-s(j)) * Q1 + s(j) * f_median; 
    Fs2(:,198+j) = (1-s(j)) * f_median + s(j) * Q3; 
    Fs2(:,297+j) = (1-s(j)) * Q3 + s(j) * maxx; 
end
d1=sqrt(trapz(t,(q_median-Q1_q).^2));
dl=sqrt(trapz(t,(Q1_q-min_q).^2));
d3=sqrt(trapz(t,(q_median-Q3_q).^2));
du=sqrt(trapz(t,(Q3_q-max_q).^2));
dmax = max(dl+d1, d3+du);
dlp = (dl+d1) / dmax;
d1p = d1 / dmax;
d3p = d3 / dmax;
dup = (d3+du) / dmax;
part1=linspace(-d1-dl,-d1,100);
part2=linspace(-d1,0,100);
part3=linspace(0,d3,100);
part4=linspace(d3,d3+du,100);
allparts=[part1,part2(2:100),part3(2:100),part4(2:100)];
[U,V]=meshgrid(t,allparts);
U=U';
V=V';
figure(316); clf;
    surf(U,V,Fs2);
    hold on;
    shading flat;
    light;
    plot3(t,zeros(1,M),f_median,'k','LineWidth',3)
    plot3(t,repmat(-d1,M,1),Q1,'b','LineWidth',3)
    plot3(t,repmat(-d1-dl,M,1),minn,'r','LineWidth',3)
    plot3(t,repmat(d3,M,1),Q3,'g','LineWidth',3)
    plot3(t,repmat(d3+du,M,1),maxx,'magenta','LineWidth',3)
    set(gca,'FontSize',18);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    
for i = 1:length(outlier_index)   % display amplitude outliers
    figure(320 + i); clf;
        for j = 1:N
            plot(t,f_tilde(:,j),'green','linewidth',2);
            hold on;
        end
        plot(t,f_median,'black','linewidth',2);
        plot(t,f_tilde(:,outlier_index(i)),'red','linewidth',2);
        xlim([t(1) t(end)]);
        ylim auto;
        set(gca,'FontSize',19)
        ti = get(gca,'TightInset');
        set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
end