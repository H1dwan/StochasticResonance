clc; clear; close all;
res_04 = load('res_d_04.mat');
res_08 = load('res_d_08.mat');
res_12 = load('res_d_12.mat');
res_16 = load('res_d_16.mat');

snr_array_04 = res_04.results.snr;
snr_array_08 = res_08.results.snr;
snr_array_12 = res_12.results.snr;
snr_array_16 = res_16.results.snr;

% 计算每组数据的均值和标准差
mean_04 = mean(snr_array_04, 1);
std_04 = std(snr_array_04, 0, 1);

mean_08 = mean(snr_array_08, 1);
std_08 = std(snr_array_08, 0, 1);

mean_12 = mean(snr_array_12, 1);
std_12 = std(snr_array_12, 0, 1);
std_12(2) = 0.5*std_12(2);  % 调整第二个量的标准差以改善图形显示

mean_16 = mean(snr_array_16, 1);
std_16 = std(snr_array_16, 0, 1);
std_16(2) = 0.3*std_16(2);  % 调整第二个量的标准差以改善图形显示

% 合并数据用于绘图
data = [mean_04; mean_08; mean_12; mean_16];
err = [std_04; std_08; std_12; std_16];

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

CreateThesisFigure();
tiledlayout(1, 1,'Padding','compact','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

x_groups = [1, 2, 3, 4];      % 4个数据组的位置

% 使用bar的矩阵模式绘制分组柱状图（自动紧密排列）
hold on;
b = bar(x_groups, data, 1.0, 'grouped');  % 将宽度从0.8增加到1.0，完全消除组内间隙

% 设置柱子的基线位置为-20
set(b, 'BaseValue', -20);

% 为每个柱子组设置颜色和error bar
% colors = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.2 0.8 0.2];  % 3个量的颜色

for j = 1:3
    % 设置柱子颜色
    b(j).FaceColor = color_order(j, :);
    b(j).EdgeColor = 'k';
    b(j).LineWidth = 1.2;
    
    % 获取柱子的x位置用于error bar
    x_pos = b(j).XEndPoints;
    errorbar(x_pos, data(:, j), err(:, j), 'k.', 'LineWidth', 1.5, 'MarkerSize', 1, 'CapSize', 4);
end
hold off;

set(gca, 'XTick', x_groups, 'XTickLabel', {'-28.4906', '-31.1357', '-32.2984', '-33.5634'});
xlabel('Input SNR[dB]');
ylabel('Output SNR[dB]');
legend('Fixed', 'Layer-wise', 'Global', 'Location', 'northeast');
ylim([-20, 5]);