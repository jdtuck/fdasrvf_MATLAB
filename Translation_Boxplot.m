function out = Translation_Boxplot(f, t)
% TRANSLATION_BOXPLOT Translation Boxplot
% -------------------------------------------------------------------------
% This function constructs the translation boxplot
%
% Usage:  out = Translation_Boxplot(f, t)
%
% Input:
% f (M,N): matrix defining N functions of M samples
% t : time vector of length M
%
% Output: structure containing
% median_t: median translation
% Q1: First quartile
% Q3: Second quartile
% minn: minimum extreme function
% maxx: maximum extreme function
% outlier_index: indexes of outlier functions


figure(100);clf;
    plot(t,f,'linewidth',2);
    set(gca,'FontSize',19);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

[~, N] = size(f);
translation = zeros(1,N);
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
    outlier_index = zeros(1,length(outlier));
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

out.time = t;
out.f = f;
out.median_t = median_t;
out.Q1 = Q1;
out.Q3 = Q3;
out.minn = minn;
out.maxx = maxx;
out.outlier_index = outlier_index;
