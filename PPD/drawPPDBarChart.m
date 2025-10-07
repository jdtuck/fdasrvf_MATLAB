function drawPPDBarChart(IndicatorMatrix, Heights, lam, idx_opt)

lam_diff = lam(2)-lam(1);
[len_lam, labelMax] = size(IndicatorMatrix);

figure

hold on
for i = 1:len_lam
    label_all_peaks = find(~isnan(Heights(i,:)));
    label_persistent_peaks = find(~isnan(IndicatorMatrix(i,:)));

    if sum(label_all_peaks) == 0
        continue
    end

    for j = 1:length(label_all_peaks)
        x = lam(i);
        y = label_all_peaks(j);

        x1 = x - lam_diff/2;
        x2 = x + lam_diff/2;
        y1 = y - 0.5;
        y2 = y + 0.5;

        if ismember(y, label_persistent_peaks)
            patch([x1, x2, x2, x1], [y1, y1, y2, y2], [0,0,0], 'EdgeColor', 'none')
        else
            patch([x1, x2, x2, x1], [y1, y1, y2, y2], [0.7,0.7,0.7], 'EdgeColor', 'none')
        end
    end
end

for j = 1:labelMax
    yline(j+0.5,'-','LineWidth',0.5)
end

xline(lam(idx_opt),'m--','linewidth',2)
hold off

xlim([lam(1),lam(end)])
ylim([0.5,labelMax+0.5])

yticks(1:labelMax);

xlabel('$\lambda$','Interpreter','latex')
ylabel('Peak Index')

set(gca,'FontSize',14, 'YDir', 'reverse');
end