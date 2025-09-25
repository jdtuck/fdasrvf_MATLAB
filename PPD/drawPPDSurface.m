function drawPPDSurface(t,lam,FNm,Heights,Locs,IndicatorMatrix,Labels,idx_opt)

n_lams = length(lam);
labelMax = size(IndicatorMatrix,2);

LocationMatrix_full = nan(n_lams,labelMax);
for i = 1:n_lams
    LocationMatrix_full(i,Labels{i}) = Locs{i};
end
LocationMatrix_sig = LocationMatrix_full .* IndicatorMatrix;
HeightMatrix_full = Heights;
HeightMatrix_sig = HeightMatrix_full .* IndicatorMatrix;

figure
surf(t,lam,FNm','LineStyle','none','FaceAlpha',0.9,'AlphaData',1,'LineWidth',0.5);
colormap(parula);

hold on
xlim([t(1),t(end)])

plot3(t,lam(idx_opt)*ones(1,length(t)),FNm(:,idx_opt),'m--','LineWidth',2)
% Loop through each label
for j = 1:labelMax

    % Find non-NaN indices for full location matrix and plot
    idx_full = find(~isnan(LocationMatrix_full(:, j)));
    plot3(t(LocationMatrix_full(idx_full, j)), ...
        lam(idx_full), ...
        HeightMatrix_full(idx_full, j), ...
        'o-', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerEdgeColor', [0.4, 0.4, 0.4], ...
        'MarkerFaceColor', [0.4, 0.4, 0.4]);

    x_offset = abs(mean(t(1) + t(LocationMatrix_full(idx_full, j)) * (t(end) - t(1))) * 0.05);
    y_offset = abs(mean(lam(idx_full)) * 0.05);
    z_offset = abs(mean(HeightMatrix_full(idx_full, j)) * 0.15);
    text(t(LocationMatrix_full(idx_full(end), j)) + x_offset, ...
        lam(idx_full(end)) + y_offset, ...
        HeightMatrix_full(idx_full(end), j) + z_offset, num2str(j), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight','bold','Color','k');

    % Find non-NaN indices for significant location matrix and plot
    idx_sig = find(~isnan(LocationMatrix_sig(:, j)));
    plot3(t(LocationMatrix_sig(idx_sig, j)), ...
        lam(idx_sig), ...
        HeightMatrix_sig(idx_sig, j), ...
        'o', 'MarkerSize', 6, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');  % Yellow markers
end

% Set axis labels with enhanced font size
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\lambda$', 'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 14);
zlabel('$\hat{g}_\lambda(t)$', 'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 14);

colormap(parula);

grid on;