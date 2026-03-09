clc; clear;close all;
%% 1. 数据导入 ===============================================
load('res_freq_stability_2.mat');

snr_005 = results(1).snr_mean;
snr_010 = results(2).snr_mean;
snr_015 = results(3).snr_mean;

mpe_005 = results(3).mpe_mean;
mpe_010 = results(2).mpe_mean;
mpe_015 = results(1).mpe_mean;

ami_005 = results(1).ami_mean;
ami_010 = results(2).ami_mean;
ami_015 = results(3).ami_mean;

rscm_005 = (1-mpe_005) .* ami_005;
rscm_010 = (1-mpe_010) .* ami_010;
rscm_015 = (1-mpe_015) .* ami_015;

D_list = results(1).D_axis;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

%% SNR 绘制 ==============================================

CreateThesisFigure(6, 4.5);
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

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mpe_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, mpe_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, mpe_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{\mathrm{MWPE}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'southeast', 'FontSize', 12);

%% AMI 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, ami_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, ami_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, ami_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{\mathrm{AMI}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);

%% MSCM 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, rscm_005, '-o', 'DisplayName', '$f=0.005$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list, rscm_010, '-s', 'DisplayName', '$f=0.010$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list, rscm_015, '-^', 'DisplayName', '$f=0.015$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{\mathrm{MSCM}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);
