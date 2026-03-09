clc; clear; close all;

SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');
ax = nexttile;
ax.Position = [0.01 0.01 0.98 0.98];

x = -1.5 : 0.01 : 1.5;
a = 1; b = 1;
plot(x, CBSR_Potential(x, a, b), 'LineWidth', 2, 'Color', color_order(1, :));
xlabel('$x$');
ylabel('$U(x)$');
% yticks([-0.2 -0.1 0 0.1]);