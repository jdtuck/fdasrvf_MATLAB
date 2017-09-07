function [median_t, Q1, Q3, minn, maxx, outlier_index] = Translation_Boxplot(f, t)

figure(100);clf;
    plot(t,f,'linewidth',2);
    set(gca,'FontSize',19);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

[M, N] = size(f);
for i = 1:N
    translation(i) = trapz(t,f(:,i)) / (t(end) - t(1));
end

figure(200);clf;
    h=boxplot(translation);
    set(h,'LineWidth',2);
    set(gca,'FontSize',19);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

h = findobj(gca,'Type','line');
m = get(h(2), 'YData');
LW = get(h(6), 'YData');
UW = get(h(7), 'YData');
minn = LW(1); Q1 = LW(2); median_t = m(1); Q3 = UW(1); maxx = UW(2);
outlier = get(h(1), 'YData');

% outliers
if isnan(outlier)
    outlier = [];
    outlier_index = [];
else
    for i = 1:length(outlier)
        outlier_index(i) = find(translation == outlier(i));
        
        figure(200+i);clf;
            plot(t,f,'green','linewidth',2);
            hold on;
            plot(t,f(:,outlier_index(i)),'red','linewidth',2);
            xlim([t(1) t(end)]);
            ylim auto;
            set(gca,'FontSize',19)
            ti = get(gca,'TightInset');
            set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    end
end