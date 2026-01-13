% =========================================================================
% Description: HSUBSR势函数的 shape factor(k1/k2) 对比以及与 CBSR、UBSR 势函数的对比
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

%% 1. 势函数参数设置
x_plot = -3:0.01:3;
a_cbsr = 1;
b_cbsr = 1;
a_ubsr = a_cbsr;
b_ubsr = b_cbsr;
xm = sqrt(a_cbsr/b_cbsr);
dU = a_cbsr^2 / (4*b_cbsr);

%% 2. 对比 CBSR、UBSR 和 HSUBSR 势函数
fig1 = CreateThesisFigure(); hold on;
y_cbsr = CBSR_Potential(x_plot, a_cbsr, b_cbsr);
y_ubsr = UBSR_Potential(x_plot, a_ubsr, b_ubsr);

shape = 2;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);

plot(x_plot, y_cbsr, '-', 'LineWidth', 2, 'DisplayName', 'CBSR Potential');
plot(x_plot, y_ubsr, '--', 'LineWidth', 2, 'DisplayName', 'UBSR Potential');
plot(x_plot, y_hsubsr, '-.', 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential ($S=2$)');
% xlim([-3 3]);
ylim([-0.25 0.25]);
xlabel('$x$');
ylabel('$U(x)$');
legend('Location', 'north');

%% 3. 不同 shape factor 下的 HSUBSR 势函数对比
fig2 = CreateThesisFigure(); hold on;
shape = 2;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);

plot(x_plot, y_hsubsr, '-', 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential ($S=2$)');

shape = 4;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, '--', 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential ($S=4$)');

shape = 10;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, '-.', 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential ($S=10$)');

shape = 50;
[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
plot(x_plot, y_hsubsr, ':', 'LineWidth', 2, 'DisplayName', 'HSUBSR Potential ($S=50$)');
xlim([-3 3]);
ylim([-0.25 0.25]);
% yticks(-0.2:0.1:0.2);
xlabel('$x$');
ylabel('$U(x)$');
legend('Location', 'north');