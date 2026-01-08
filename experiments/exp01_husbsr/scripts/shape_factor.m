% =========================================================================
% Description: HSUBSR势函数的 shape factor(k1/k2) 对比
%
% Author: LiuShuang
% Created: 2025-12-15
% Last Modified: 2025-12-15
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

SetThesisDefaultStyle();
fig = CreateThesisFigure();
a_ubsr = 1;
b_ubsr = 1;
xm = sqrt(a_ubsr/b_ubsr);
dU = a_ubsr^2 / (4*b_ubsr);

x_plot = -3:0.01:3;
y_ubsr = UBSR_Potential(x_plot, 'a', a_ubsr, 'b', b_ubsr);

shape = 2;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);

plot(x_plot, y_ubsr, 'LineWidth', 2, 'DisplayName', 'UBSR Potential'); hold on
plot(x_plot, y_hsubsr, 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential (shape=2)');

shape = 4;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential (shape=4)');

shape = 6;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential (shape=6)');

shape = 50;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential (shape=50)');

ylim([-0.25 0.25]);
xlabel('x');
ylabel('U(x)');
legend('Location', 'best');