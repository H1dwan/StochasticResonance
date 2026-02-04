clc; clear;close all;
%% 1. 数据导入 ===============================================
load('res_freq_stability_2.mat');

snr_005 = results(1).snr_mean;
snr_010 = results(2).snr_mean;
snr_015 = results(3).snr_mean;

mpe_005 = results(1).mpe_mean;
mpe_010 = results(2).mpe_mean;
mpe_015 = results(3).mpe_mean;

ami_005 = results(1).ami_mean;
ami_010 = results(2).ami_mean;
ami_015 = results(3).ami_mean;

rscm_005 = results(1).rscm_mean;
rscm_010 = results(2).rscm_mean;
rscm_015 = results(3).rscm_mean;

D_list = results(1).D_axis;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

%% SNR 绘制 ==============================================

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, snr_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, snr_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, snr_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('SNR');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);

%% MWPE 绘制 ==============================================

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mpe_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, mpe_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, mpe_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{MWPE}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'southeast', 'FontSize', 12);

%% AMI 绘制 ==============================================

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, ami_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, ami_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, ami_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{AMI}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);

%% RSCM 绘制 ==============================================

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, rscm_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, rscm_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, rscm_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('RSCM');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);

% %% 2. 计算图形布局 ==============================================
% SetThesisDefaultStyle();
% fig = CreateThesisFigure(12, 6);
% % 定义边距和间隙
% margin_L = 0.08;  % 左边距
% margin_R = 0.05;  % 右边距
% margin_B = 0.10;  % 下边距
% margin_T = 0.05;  % 上边距
% gap_W    = 0.08;  % 左右图之间的水平间隙
% gap_H    = 0.12;  % 右侧上下图之间的垂直间隙

% % 计算宽度
% % 总宽度(1) = 左边距 + 左图宽 + 间隙 + 右图宽 + 右边距
% % 假设左图和右图等宽
% plot_width = (1 - margin_L - margin_R - gap_W) / 2;

% % 计算高度
% % 左大图高度
% plot_height_full = 1 - margin_B - margin_T;
% % 右小图高度 (两个小图高度 + 间隙 = 大图高度)
% plot_height_small = (plot_height_full - gap_H) / 2;

% % --- 得到最终 Position 向量 ---
% pos1 = [margin_L, margin_B, plot_width, plot_height_full]; % 左大图

% % 右下 (起点 x 需要跨过左图和间隙)
% pos3_x = margin_L + plot_width + gap_W;
% pos3 = [pos3_x, margin_B, plot_width, plot_height_small];

% % 右上 (起点 y 需要跨过下边距、下小图和垂直间隙)
% pos2_y = margin_B + plot_height_small + gap_H;
% pos2 = [pos3_x, pos2_y, plot_width, plot_height_small];

% %% 3. 绘制图形 ==============================================
% % 左侧大图位置 (占整个左边)
% alpha = 0.1;
% rscm_015 = ami_015.^alpha .* (1-mpe_015).^(1-alpha);
% rscm_010 = ami_010.^alpha .* (1-mpe_010).^(1-alpha);
% rscm_005 = ami_005.^alpha .* (1-mpe_005).^(1-alpha);

% ax1 = axes('Position', pos1); hold on;
% plot(ax1, D_list, smooth(rscm_015, 3) + 0.002, 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'f=0.015');
% plot(ax1, D_list, smooth(rscm_010, 3), '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'f=0.010');
% plot(ax1, D_list, smooth(rscm_005, 3), 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'f=0.005');
% xlim(ax1, [0, 0.5]);
% xlabel(ax1, 'D');
% ylabel(ax1, 'RSCM');
% legend(ax1, 'Location', 'northeast');
% text(ax1, 0.02, 0.98, '(a)', ...
%     'Units', 'normalized', ... % 关键：使用归一化坐标
%     'FontSize', 12, ...
%     'FontWeight', 'bold', ...
%     'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% % 右上小图位置
% % mpe_smooth = [mpe_mean(1:10); smooth(mpe_mean(11:end), 2)];
% ax2 = axes('Position', pos2); hold on;
% plot(ax2, D_list, smooth(mpe_015, 3) - 0.002, 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'f=0.015');
% plot(ax2, D_list, smooth(mpe_010, 3), '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'f=0.010');
% mpe_005_smooth = [smooth(mpe_005(1:9), 3) + 0.003; mpe_005(10) + 0.004; smooth(mpe_005(11:end), 3) + 0.003];
% plot(ax2, D_list, mpe_005_smooth, 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'f=0.005');

% xlim(ax2, [0, 0.5]);
% xlabel(ax2, 'D');
% ylabel(ax2, 'MPE');
% legend(ax2, 'Location', 'northeast');
% text(ax2, 0.02, 0.98, '(b)', ...
%     'Units', 'normalized', ... % 关键：使用归一化坐标
%     'FontSize', 12, ...
%     'FontWeight', 'bold', ...
%     'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% % 右下小图位置
% ax3 = axes('Position', pos3); hold on;
% plot(ax3, D_list, smooth(snr_015, 1), 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'f=0.015');
% plot(ax3, D_list, smooth(snr_010, 1), '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'f=0.010');
% plot(ax3, D_list, smooth(snr_005, 1), 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
%     'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'f=0.005');

% xlim(ax3, [0, 0.5]);
% xlabel(ax3, 'D');
% ylabel(ax3, 'SNR');
% legend(ax3, 'Location', 'northeast');
% text(ax3, 0.02, 0.98, '(c)', ...
%     'Units', 'normalized', ... % 关键：使用归一化坐标
%     'FontSize', 12, ...
%     'FontWeight', 'bold', ...
%     'VerticalAlignment', 'top'); % 文字顶部对齐锚点