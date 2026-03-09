clc; clear; close all;

load('res_ami_vs_d.mat');
snr_list = results.snr;
ami_list = results.ami;
mwpe_list = results.mpe;
D_list = results.D_list;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

%% 绘制 SNR 和 AMI 的双坐标图 ==============================================
CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];
yyaxis left;
plot(D_list, snr_list, '-o', 'DisplayName', 'SNR', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
ylabel('SNR');
yyaxis right;
plot(D_list, ami_list, '-s', 'DisplayName', '$J_{\mathrm{AMI}}$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
ylabel('$J_{\mathrm{AMI}}$');

xlim([0.0, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);


%% 绘制 SNR 和 MWPE 的双坐标图 ==============================================
CreateThesisFigure(6, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];
yyaxis left;
plot(D_list, snr_list, '-o', 'DisplayName', 'SNR', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
ylabel('SNR');
yyaxis right;
plot(D_list, mwpe_list, '-s', 'DisplayName', '$J_{\mathrm{MWPE}}$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
ylabel('$J_{\mathrm{MWPE}}$');

xlim([0.0, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'southeast', 'FontSize', 12);