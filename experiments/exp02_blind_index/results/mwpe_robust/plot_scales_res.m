clc;clear; close all;

%% 有效性验证实验
% load('res_m_3_d_014_t_mean_r_857_924.mat');
% snr_list = results.output.snr_list;
% mpe_list = results.output.mpe_list;
% D_list = 0.05:0.01:0.45;

% SetThesisDefaultStyle();
% CreateThesisFigure();
% tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

% color_order = get(0, 'DefaultAxesColorOrder');

% ax = nexttile;
% ax.Position = [0.01 0.01 0.98 0.98];
% yyaxis left;
% plot(D_list, snr_list, '-o', 'DisplayName', 'SNR', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
% ylabel('SNR');
% yyaxis right;
% plot(D_list, mpe_list, '-s', 'DisplayName', '$J_{MWPE}$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
% ylabel('$J_{MWPE}$');


%% 嵌入维数 m 稳健性实验
% load('res_m_345_t_mean_2.mat');

% mpe_list_3 = results.output.mpe_list;
% mpe_list_4 = results.output.mpe_list4;
% mpe_list_5 = results.output.mpe_list5;
% snr_list = results.output.snr_list;
% D_list = 0.05:0.01:0.45;

% SetThesisDefaultStyle();
% CreateThesisFigure();
% tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

% color_order = get(0, 'DefaultAxesColorOrder');

% ax = nexttile;
% ax.Position = [0.01 0.01 0.98 0.98];

% plot(D_list, mpe_list_3, '-o', 'DisplayName', '$m=3$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
% hold on;
% plot(D_list, mpe_list_4, '-s', 'DisplayName', '$m=4$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
% plot(D_list, mpe_list_5, '-^', 'DisplayName', '$m=5$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
% ylabel('$J_{MWPE}$');
% xlabel('$D$');

% xlim([0.00, 0.50]);
% xticks(0:0.1:0.50);
% legend('Location', 'northeast', 'FontSize', 12);

% CreateThesisFigure();
% tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');
% ax = nexttile;
% ax.Position = [0.01 0.01 0.98 0.98];
% plot(D_list, snr_list, '-o', 'DisplayName', 'SNR', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
% ylabel('SNR');
% xlabel('$D$');

% xlim([0.00, 0.50]);
% xticks(0:0.1:0.50);
% legend('Location', 'northeast', 'FontSize', 12);

%% 映射系数 C 稳健性实验
% load('res_t_mean_c_246.mat');
% mpe_list_2 = results.output.mpe_list2;
% mpe_list_4 = results.output.mpe_list4;
% mpe_list_6 = results.output.mpe_list6;
% load('res_t_mean_c_50.mat');
% mpe_list_50 = results.output.mpe_list2;
% D_list = 0.05:0.01:0.45;


% SetThesisDefaultStyle();
% CreateThesisFigure();
% tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

% color_order = get(0, 'DefaultAxesColorOrder');

% ax = nexttile;
% ax.Position = [0.01 0.01 0.98 0.98];

% plot(D_list, mpe_list_2, '-o', 'DisplayName', '$C=2$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
% hold on;
% plot(D_list, mpe_list_4, '-s', 'DisplayName', '$C=4$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
% plot(D_list, mpe_list_6, '-^', 'DisplayName', '$C=6$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
% plot(D_list, mpe_list_50, '-d', 'DisplayName', '$C=50$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(4, :));

% ylabel('$J_{MWPE}$');
% xlabel('$D$');

% xlim([0.00, 0.50]);
% xticks(0:0.1:0.50);
% legend('Location', 'northeast', 'FontSize', 12);

%% 融合策略对比实验
load('res_m_345_t_mean_2.mat');

mpe_list_3_mean = results.output.mpe_list;
mpe_list_4_mean = results.output.mpe_list4;
mpe_list_5_mean = results.output.mpe_list5;

load('res_t_min_m_345.mat')
mpe_list_3_min = results.output.mpe_list3;
mpe_list_4_min = results.output.mpe_list4;
mpe_list_5_min = results.output.mpe_list5;

D_list = 0.05:0.01:0.45;

SetThesisDefaultStyle();
CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

color_order = get(0, 'DefaultAxesColorOrder');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mpe_list_3_mean, '-o', 'DisplayName', '$m=3$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, mpe_list_4_mean, '-s', 'DisplayName', '$m=4$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, mpe_list_5_mean, '-^', 'DisplayName', '$m=5$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{MWPE}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

color_order = get(0, 'DefaultAxesColorOrder');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mpe_list_3_min, '-o', 'DisplayName', '$m=3$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, mpe_list_4_min, '-s', 'DisplayName', '$m=4$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, mpe_list_5_min, '-^', 'DisplayName', '$m=5$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{MWPE}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);