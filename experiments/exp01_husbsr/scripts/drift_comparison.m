% =========================================================================
% Description: 比较三种双稳态势函数及其漂移项 CBSR UBSR HSUBSR
%
% Author: LiuShuang
% Created: 2025-12-07
% Last Modified: 2025-12-18
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

%% 1. 势函数设置 =============================================
% CBSR 势函数参数
a_cbsr = 1;
b_cbsr = 1;
xm = sqrt(a_cbsr / b_cbsr);
dU = a_cbsr^2 / (4 * b_cbsr);
drift_cbsr = @(x) CBSR_Dynamics(x, a_cbsr, b_cbsr);

% UBSR 势函数参数
a_ubsr = 1;
b_ubsr = 1;
drift_ubsr = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);

% HSUBSR 势函数参数
shape_factor = 2; % shape factor k1/k2
[a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr] = CalibrateHSUBSR(xm, dU, shape_factor);
drift_hsubsr = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);

%% 2. 势函数及漂移项对比 =============================================
fig = CreateThesisFigure(14, 6);
subplot(1,2,1);
x = -3:0.01:3;
plot(x, CBSR_Potential(x, a_cbsr, b_cbsr), '-', 'LineWidth', 2); hold on;
plot(x, UBSR_Potential(x, a_ubsr, b_ubsr), '--', 'LineWidth', 2);
plot(x, HSUBSR_Potential(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr)-...
    HSUBSR_Potential(0, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr), '-.', 'LineWidth', 2);
ylim([-dU, dU]);
xlabel('$x$')
ylabel('$U(x)$')
legend('CBSR', 'UBSR', 'HSUBSR', 'Location', 'north');

subplot(1,2,2);
plot(x, drift_cbsr(x), '-', 'LineWidth', 2); hold on;
plot(x, drift_ubsr(x), '--', 'LineWidth', 2);
plot(x, drift_hsubsr(x), '-.', 'LineWidth', 2);
legend('CBSR', 'UBSR', 'HSUBSR', 'Location', 'northeast');
ylim([-1, 1]);
xlabel('$x$')
ylabel('-$U^{\prime}(x)$')