clc; clear; close all;
%% 1. 数据导入 ===============================================
load('res_amp_stability.mat');

snr_005 = results(1).snr_mean;
snr_010 = results(2).snr_mean;
snr_015 = results(3).snr_mean;

mpe_005 = results(1).mpe_mean;
mpe_010 = results(2).mpe_mean;
mpe_015 = results(3).mpe_mean;

ami_005 = results(1).ami_mean;
ami_010 = results(2).ami_mean;
ami_015 = results(3).ami_mean;

mscm_005 = results(1).rscm_mean;
mscm_010 = results(2).rscm_mean;
mscm_015 = results(3).rscm_mean;

D_list = 0.05:0.01:0.45;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

%% SNR 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list(3:end), snr_005(3:end), '-o', 'DisplayName', '$A=0.05$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list(3:end), snr_010(3:end), '-s', 'DisplayName', '$A=0.10$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list(3:end), snr_015(3:end), '-^', 'DisplayName', '$A=0.15$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
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

plot(D_list(3:end), mpe_005(3:end), '-o', 'DisplayName', '$A=0.05$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list(3:end), mpe_010(3:end), '-s', 'DisplayName', '$A=0.10$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list(3:end), mpe_015(3:end), '-^', 'DisplayName', '$A=0.15$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
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

plot(D_list(3:end), ami_005(3:end), '-o', 'DisplayName', '$A=0.05$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list(3:end), ami_010(3:end), '-s', 'DisplayName', '$A=0.10$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list(3:end), ami_015(3:end), '-^', 'DisplayName', '$A=0.15$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
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

plot(D_list(3:end), mscm_005(3:end), '-o', 'DisplayName', '$A=0.05$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
hold on;
plot(D_list(3:end), mscm_010(3:end), '-s', 'DisplayName', '$A=0.10$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
plot(D_list(3:end), mscm_015(3:end), '-^', 'DisplayName', '$A=0.15$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{\mathrm{MSCM}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);
