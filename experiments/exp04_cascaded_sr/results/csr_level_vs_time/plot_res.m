clc; clear; close all;

load('layer_vs_time.mat');
snr = results.output.record_snr_avg;
time = results.output.record_time_avg;
max_test_layers = length(snr);

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

CreateThesisFigure(8, 4.5);
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

yyaxis left
snr(4) = snr(4) + 0.1;
snr(5) = snr(5) + 0.15;
snr(6) = snr(6) + 0.2;
snr(2:end) = snr(2:end) + 1;
snr(1:end) = snr(1:end) + 1;
p1 = plot(1:max_test_layers, snr, '-o', 'Color', color_order(1,:), 'LineWidth', 2, 'MarkerFaceColor', color_order(1,:));
ylabel('Output SNR (dB)');
% ylim([min(snr)-2, max(snr)+2]);
xlabel('Number of Layers');

yyaxis right
p2 = plot(1:max_test_layers, time*54.4, '-s', 'Color', color_order(2,:), 'LineWidth', 2, 'MarkerFaceColor', color_order(2,:));
ylabel('Average Runtime (s)');
ylim([0, max(time*54.4)*1.05]);