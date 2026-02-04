clc; clear; close all;

load('res_ami_vs_d_2.mat');
snr_list = results.output.snr_list;
ami_list = results.output.ami_list;
D_list = 0.07:0.01:0.45;

SetThesisDefaultStyle();
CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

color_order = get(0, 'DefaultAxesColorOrder');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];
yyaxis left;
plot(D_list, snr_list, '-o', 'DisplayName', 'SNR', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(1, :));
ylabel('SNR');
yyaxis right;
plot(D_list, ami_list, '-s', 'DisplayName', '$J_{AMI}$', 'LineWidth', 1.5, 'MarkerFaceColor', color_order(2, :));
ylabel('$J_{AMI}$');

xlim([0.0, 0.50]);
xticks(0:0.1:0.50);
legend('Location', 'northeast', 'FontSize', 12);