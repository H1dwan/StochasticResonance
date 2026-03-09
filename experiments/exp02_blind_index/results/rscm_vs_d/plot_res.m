clc; clear; close all;

load('res_rscm_vs_d .mat');
snr_list = results.snr;
ami_list = results.ami;
mwpe_list = results.mpe;
mscm_list = results.rscm;
D_list = results.D_list;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

%% SNR 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, snr_list, '-o', 'LineWidth', 1.5, 'Color', color_order(1, :), 'MarkerFaceColor', color_order(1, :));
ylabel('SNR');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);

%% MWPE 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mwpe_list, '-^', 'LineWidth', 1.5, 'Color', color_order(2, :), 'MarkerFaceColor', color_order(2, :));
ylabel('$J_{\mathrm{MWPE}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);

%% AMI 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, ami_list, '-s', 'LineWidth', 1.5, 'Color', color_order(3, :), 'MarkerFaceColor', color_order(3, :));
ylabel('$J_{\mathrm{AMI}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);

%% MSCM 绘制 ==============================================

CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, mscm_list, '-d',  'LineWidth', 1.5, 'Color', color_order(4, :), 'MarkerFaceColor', color_order(4, :));
ylabel('$J_{\mathrm{MSCM}}$');
xlabel('$D$');

xlim([0.00, 0.50]);
xticks(0:0.1:0.50);