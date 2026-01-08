clc;clear;
%% 1. 数据导入 ===============================================
D_list = 0.05:0.01:0.45;

results_015 = load('results_015amp.mat');
% snr_015 = results_015.results.snr_mean;
[snr_015, ~, ~] = HSUBSR_SNR(1, 0.25, 50, 0.010, 0.15, D_list);
mpe_015 = results_015.results.mpe_mean;
ami_015 = results_015.results.ami_mean;

results_010 = load('results_010amp.mat');
% snr_010 = results_010.results.snr_mean;
[snr_010, ~, ~] = HSUBSR_SNR(1, 0.25, 50, 0.010, 0.10, D_list);
mpe_010 = results_010.results.mpe_mean;
ami_010 = results_010.results.ami_mean;

reuslts_005 = load('results_005amp.mat');
% snr_005 = reuslts_005.results.snr_mean;
[snr_005, ~, ~] = HSUBSR_SNR(1, 0.25, 50, 0.010, 0.05, D_list);
mpe_005 = reuslts_005.results.mpe_mean;
ami_005 = reuslts_005.results.ami_mean;



%% 2. 计算图形布局 ==============================================
SetThesisDefaultStyle();
fig = CreateThesisFigure(12, 6);
% 定义边距和间隙
margin_L = 0.08;  % 左边距
margin_R = 0.05;  % 右边距
margin_B = 0.10;  % 下边距
margin_T = 0.05;  % 上边距
gap_W    = 0.08;  % 左右图之间的水平间隙
gap_H    = 0.12;  % 右侧上下图之间的垂直间隙

% 计算宽度
% 总宽度(1) = 左边距 + 左图宽 + 间隙 + 右图宽 + 右边距
% 假设左图和右图等宽
plot_width = (1 - margin_L - margin_R - gap_W) / 2;

% 计算高度
% 左大图高度
plot_height_full = 1 - margin_B - margin_T;
% 右小图高度 (两个小图高度 + 间隙 = 大图高度)
plot_height_small = (plot_height_full - gap_H) / 2;

% --- 得到最终 Position 向量 ---
pos1 = [margin_L, margin_B, plot_width, plot_height_full]; % 左大图

% 右下 (起点 x 需要跨过左图和间隙)
pos3_x = margin_L + plot_width + gap_W;
pos3 = [pos3_x, margin_B, plot_width, plot_height_small];

% 右上 (起点 y 需要跨过下边距、下小图和垂直间隙)
pos2_y = margin_B + plot_height_small + gap_H;
pos2 = [pos3_x, pos2_y, plot_width, plot_height_small];

%% 3. 绘制图形 ==============================================
% 左侧大图位置 (占整个左边)
alpha = 0.1;
rscm_015 = ami_015.^alpha .* (1-mpe_015).^(1-alpha);
rscm_010 = ami_010.^alpha .* (1-mpe_010).^(1-alpha);
rscm_005 = ami_005.^alpha .* (1-mpe_005).^(1-alpha);

ax1 = axes('Position', pos1); hold on;
plot(ax1, D_list, smooth(rscm_015, 3) + 0.005, 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'A=0.15');
plot(ax1, D_list, smooth(rscm_010, 3), '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'A=0.10');
plot(ax1, D_list, smooth(rscm_005, 3) - 0.005, 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'A=0.05');
xlim(ax1, [0, 0.5]);
xlabel(ax1, 'D');
ylabel(ax1, 'RSCM');
legend(ax1, 'Location', 'northeast');
text(ax1, 0.02, 0.98, '(a)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% 右上小图位置
ax2 = axes('Position', pos2); hold on;
mpe_015_smooth = [mpe_015(1:12); smooth(mpe_015(13:end), 4)];
plot(ax2, D_list, mpe_015_smooth - 0.005, 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'A=0.15');

mpe_010_smooth = [mpe_010(1:11); smooth(mpe_010(12:end), 4)];
plot(ax2, D_list, mpe_010_smooth, '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'A=0.10');

mpe_005_smooth = [mpe_005(1:11); smooth(mpe_005(12:end), 4)];
plot(ax2, D_list, mpe_005_smooth + 0.005, 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'A=0.05');

xlim(ax2, [0, 0.5]);
xlabel(ax2, 'D');
ylabel(ax2, 'MPE');
legend(ax2, 'Location', 'north');
text(ax2, 0.02, 0.98, '(b)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% 右下小图位置
ax3 = axes('Position', pos3); hold on;
plot(ax3, D_list, smooth(snr_015, 1), 's-', 'Color', [0.0000, 0.6196, 0.4510], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'A=0.15');
plot(ax3, D_list, smooth(snr_010, 1), '^-', 'Color', [0.2510, 0.3686, 0.6353], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'A=0.10');
plot(ax3, D_list, smooth(snr_005, 1), 'o-', 'Color', [0.9020, 0.2941, 0.2078], ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'A=0.05');
xlim(ax3, [0, 0.5]);
xlabel(ax3, 'D');
ylabel(ax3, 'SNR');
legend(ax3, 'Location', 'northeast');
text(ax3, 0.02, 0.98, '(c)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点