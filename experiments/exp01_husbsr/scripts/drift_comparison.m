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

% % PLBSR 势函数参数
% a_plbsr = 1.7;
% b_plbsr = 1;
% c_plbsr = 0.25;
% drift_plbsr = @(x) PLBSR_Dynamics(x, a_plbsr, b_plbsr, c_plbsr);

% HSUBSR 势函数参数
shape_factor = 20; % shape factor k1/k2
[a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr] = CalibrateHSUBSR(xm, dU, shape_factor);
drift_hsubsr = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);

%% 2. 势函数及漂移项对比 =============================================
SetThesisDefaultStyle();
fig = CreateThesisFigure(16, 6);
subplot(1,2,1);
x = -4:0.01:4;
plot(x, CBSR_Potential(x, 'a', a_cbsr, 'b', b_cbsr), '--'); hold on;
plot(x, UBSR_Potential(x, 'a', a_ubsr, 'b', b_ubsr), ':'); 
% plot(x, PLBSR_Potential(x, 'a', a_plbsr, 'b', b_plbsr, 'c', c_plbsr), ':');
plot(x, HSUBSR_Potential(x, 'a', a_hsubsr, 'b', b_hsubsr, 'k1', k1_hsubsr, 'k2', k2_hsubsr)-...
    HSUBSR_Potential(0, 'a', a_hsubsr, 'b', b_hsubsr, 'k1', k1_hsubsr, 'k2', k2_hsubsr), '-.');
ylim([-dU, dU]);
legend('CBSR', 'UBSR', 'HSUBSR', 'Location', 'north');

subplot(1,2,2);
plot(x, drift_cbsr(x), '--'); hold on;
plot(x, drift_ubsr(x), ':');
plot(x, drift_hsubsr(x), '-.');
legend('CBSR', 'UBSR', 'HSUBSR', 'Location', 'northeast');

% ExportThesisFigure(fig, 'Potential_Drift');