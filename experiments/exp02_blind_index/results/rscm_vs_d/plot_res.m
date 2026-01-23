clc;clear;close all;

load("mpe_14.mat");
load("ami_13.mat");
% load("snr_13.mat");
% snr = snr_mean;
load("snr_theo.mat");
% snr = results.snr_theo(1:41);

D_list = 0.05:0.01:0.45;

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

% 左侧大图位置 (占整个左边)
alpha = 0.1;
% rscm = ami_mean.^alpha .* (1-mpe_mean).^(1-alpha);
rscm = ami_mean .* (1-mpe_mean-0.745);
CreateThesisFigure(8, 6);
layout = tiledlayout(1,1);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98]; 
% rscm_smooth = [rscm(1:9); smooth(rscm(10:end), 2)];
plot(D_list, rscm, 'o-', 'Color', [0.9020, 0.2941, 0.2078], 'MarkerSize', 5, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'RSCM');
xlim([0, 0.5]);
xlabel('$D$');
ylabel('RSCM');
xticks(0:0.1:0.5);
yticks(0.01:0.01:0.06);
text(0.02, 0.98, '(a)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% ax1 = axes('Position', pos1);
plot(ax1, D_list, smooth(rscm, 1), 'o-', 'Color', [0.9020, 0.2941, 0.2078], 'MarkerSize', 5, 'MarkerFaceColor', [0.9020, 0.2941, 0.2078], 'DisplayName', 'RSCM');
xlim(ax1, [0, 0.5]);
xlabel(ax1, 'D');
ylabel(ax1, 'RSCM');
text(ax1, 0.02, 0.98, '(a)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% 右上小图位置
CreateThesisFigure(8, 6);
layout = tiledlayout(1,1);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98]; 

mpe_smooth = [mpe_mean(1:10); smooth(mpe_mean(11:end), 3)];
plot(D_list, mpe_smooth + 0.745, '^-', 'Color', [0.2510, 0.3686, 0.6353], 'MarkerSize', 5, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'MPE');
xlim([0, 0.5]);
xticks(0:0.1:0.5);
yticks(0.75:0.01:0.80);
% yticklabels({'0.74', '0.75', '0.76', '0.77', '0.78', '0.79', '0.80'});
xlabel('$D$');
ylabel('MWPE');

ax2 = axes('Position', pos2);
plot(ax2, D_list, mpe_smooth, '^-', 'Color', [0.2510, 0.3686, 0.6353], 'MarkerSize', 5, 'MarkerFaceColor', [0.2510, 0.3686, 0.6353], 'DisplayName', 'MPE');
xlim(ax2, [0, 0.5]);
xlabel(ax2, 'D');
ylabel(ax2, 'MPE');
text(ax2, 0.02, 0.98, '(b)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点

% 右下小图位置
CreateThesisFigure(8, 6);
layout = tiledlayout(1,1);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98]; 
[~,snr,~] = HSUBSR_SNR(1, 0.25, 20, 0.01, 0.1, 0:0.01:0.5);
plot(0:0.01:0.5, snr, 's-', 'Color', [0.0000, 0.6196, 0.4510], 'MarkerSize', 5, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'SNR');
xlim([0, 0.5]);
xticks(0:0.1:0.5);
yticks(0:0.01:0.05);
% ylim(ax3, [0, 0.017]);
xlabel('$D$');
ylabel('SNR');


ax3 = axes('Position', pos3);
plot(ax3, D_list, smooth(snr, 3), 's-', 'Color', [0.0000, 0.6196, 0.4510], 'MarkerSize', 5, 'MarkerFaceColor', [0.0000, 0.6196, 0.4510], 'DisplayName', 'SNR');
xlim(ax3, [0, 0.5]);
% ylim(ax3, [0, 0.017]);
xlabel(ax3, '$D$');
ylabel(ax3, 'SNR');
text(ax3, 0.02, 0.98, '(c)', ...
    'Units', 'normalized', ... % 关键：使用归一化坐标
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top'); % 文字顶部对齐锚点             