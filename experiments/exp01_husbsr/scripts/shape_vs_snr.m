% =========================================================================
% Description: 研究在固定势垒高度和势间距下，shape factor 对 HSUBSR 输出 SNR
%              的影响
%              dU = 0.2 : [6 20 50]
%              dU = 0.3 : [4 20 50]
%              dU = 0.4 : [4 6 20]
%              dU = 0.5 : [4 50]
%              dU = 0.6 : [50]
%              dU = 0.7 : [4]
%
% Author: LiuShuang
% Created: 2025-12-24
% Last Modified: 2025-12-24
%
% Caution: Snr 不是随着 shape factor 的增大而单调变化的，是先增大后减小再增大的
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;
SetThesisDefaultStyle();

%% 1. 实验设置 ============================================================
% --- 锁定的几何目标 (控制变量) ---
xm = 1;
dU = 1;

% --- 待对比的 Shape Factors ---
shape_list = [1.5, 2, 6, 20];

% --- 信号与系统参数 ---
A = 0.05;                  % 弱信号幅值
f0 = 0.01;                % 信号频率
fs = 5;                   % 采样率
T_sim = 2000;             % 仿真时长
dt = 1/fs;
t = (0:dt:T_sim-dt)';
N = length(t);

% --- 噪声扫描 ---
D_list = 0.0:0.01:2;     % 关注低噪声到中噪声区
n_repeat = 20;            % 重复次数取平均

%% 2. 仿真主循环 ==========================================================
s_clean = A * sin(2*pi*f0*t); % 纯净信号
results = struct();

figure('Position', [100, 100, 1000, 400], 'Color', 'w');

% --- 子图1: 势函数形状对比 ---
subplot(1, 2, 1); hold on;
x_plot = linspace(-3, 3, 1000);

for k = 1:length(shape_list)
    shape = shape_list(k);
    [a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
    y_hsubsr = HSUBSR_Potential(x_plot, a, b, k1, k2);
    plot(x_plot, y_hsubsr, 'LineWidth', 2, 'DisplayName', sprintf('HSUBSR (shape=%d)', shape));
    
    results(k).shape = shape;
    results(k).params = [a, b, k1, k2];
end
legend('Location', 'best');

% --- 子图2：信噪比对比 ---
subplot(1, 2, 2); hold on;
linestyle_list = {'-', '--', '-.', ':'};
for k = 1:length(shape_list)
    shape = results(k).shape;
    params = results(k).params;
    [~, snr_adiabatic, ~] = HSUBSR_SNR(xm, dU, shape, f0, A, D_list);
    plot(D_list, snr_adiabatic, linestyle_list{k}, 'LineWidth', 1.5, 'DisplayName', sprintf('Theory (shape=%d)', shape));
end


% for k = 1:length(shape_list)
%     shape = results(k).shape;
%     params = results(k).params;
%     SNR_out = zeros(size(D_list));
%     drift = @(x) HSUBSR_Dynamics(x, params(1), params(2), params(3), params(4));
%     fprintf('Simulating for shape factor = %d...\n', shape);

%     for i = 1:length(D_list)
%         D = D_list(i);
%         snr_accum = 0;

%         parfor n = 1:n_repeat
%             noise = sqrt(2*D*fs) * randn(N, 1);
%             x_out = RK4Solver2(drift, s_clean, noise, fs);
%             snr_accum = snr_accum + SNRo2(x_out, fs, f0);
%         end

%         SNR_out(i) = snr_accum / n_repeat;
%     end

%     plot(D_list, SNR_out, '-o', 'LineWidth', 2, 'DisplayName', sprintf('Shape=%d', shape));
% end

legend('Location', 'best');