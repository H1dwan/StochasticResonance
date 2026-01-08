% =========================================================================
% Description: 轨迹均方误差(trajectory RMS error)对比实验
%              比较 HSUBSR / UBSR / PLBSR 在不同采样率下
%              对精细参考解的逼近误差，用于度量数值稳定性
%
% Author: LiuShuang
% Created: 2025-12-05
% Last Modified: 2025-12-22
%
% Usage:
%   使用保存的数据进行绘图
% =========================================================================

clc; clear; close all;
rng(2025);   % 固定随机种子，保证可重复性

%% 1. 全局仿真参数 =====================================================

t_total      = 800;      % 总仿真时间 (s)
f0           = 0.01;     % 弱周期信号频率 (Hz)
a0           = 0.05;     % 弱周期信号幅值
d_noise      = 0.25;     % 噪声强度 D

fs_ref       = 1000;                   % 参考采样率 (高采样，高精度解)
fs_list      = [10, 20, 25, 40, 50, 100, 125, 200];      % 需要比较的采样率 (必须整除 fs_ref)
num_fs       = numel(fs_list);

% -------- 固定的势函数结构参数 ----------------------------------------
xm = 1;
dU = 0.25;

% -------- HSUBSR 势函数参数 ------------------------------------------
shape_factor = 4;
[a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr] = CalibrateHSUBSR(xm, dU, shape_factor);

% -------- UBSR 势函数参数 --------------------------------------------
a_ubsr     = 1.0;
b_ubsr     = 1.0;

% -------- PLBSR 势函数参数 -------------------------------------------
u_plbsr    = 0.25;
l_plbsr    = 1;

%% 2. 构造参考时间栅格、信号和布朗运动路径 =============================

h_ref        = 1 / fs_ref;
n_steps_ref  = t_total * fs_ref;          % 参考步数
n_ref        = n_steps_ref + 1;           % 参考采样点数
t_ref        = (0:n_ref-1)' * h_ref;

% 参考信号 (高采样率)
clean_signal_ref = a0 * sin(2*pi*f0*t_ref);

% 生成一条布朗运动路径 dW_ref
% dW_ref(i) ~ N(0, h_ref)，累计和是标准布朗运动
dW_ref = sqrt(h_ref) * randn(n_steps_ref, 1);

% 对应的随机力序列 noise_seq_ref，使得：
%   dx_stoch = noise_seq_ref(i) * h_ref = sqrt(2D) dW_ref(i)
%   => noise_seq_ref(i) = sqrt(2D) * dW_ref(i) / h_ref
noise_seq_ref = zeros(n_ref, 1);
noise_seq_ref(1:n_steps_ref) = sqrt(2 * d_noise) * dW_ref / h_ref;
noise_seq_ref(end) = 0;    % 最后一个不会被用到，仅占位

%% 3. 计算三种模型在参考采样率下的轨迹 (近似真解) =======================

drift_hsubsr = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);
drift_ubsr = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
drift_plbsr  = @(x) PLBSR_Dynamics(x, u_plbsr, l_plbsr);

x_ref_hsubsr = RK4Solver2(drift_hsubsr, clean_signal_ref, noise_seq_ref, fs_ref);
x_ref_ubsr = RK4Solver2(drift_ubsr, clean_signal_ref, noise_seq_ref, fs_ref);
x_ref_plbsr  = RK4Solver2(drift_plbsr,  clean_signal_ref, noise_seq_ref, fs_ref);
%% 4. 不同采样率下轨迹误差计算 =========================================

rms_err_hsubsr = zeros(num_fs, 1);
rms_err_ubsr = zeros(num_fs, 1);
rms_err_plbsr  = zeros(num_fs, 1);

for i_fs = 1:num_fs
    fs_coarse      = fs_list(i_fs);
    h_coarse       = 1 / fs_coarse;
    
    % 要求 fs_ref 能被 fs_coarse 整除
    if mod(fs_ref, fs_coarse) ~= 0
        error('fs_ref 必须能被 fs_list 中的每个采样率整除！');
    end
    
    m = fs_ref / fs_coarse;              % 一个粗步长包含的细步长个数
    
    % 粗采样步数和采样点数
    n_steps_coarse = t_total * fs_coarse;
    n_coarse       = n_steps_coarse + 1;
    
    % 粗采样时间栅格与信号
    t_coarse            = (0:n_coarse-1)' * h_coarse;
    clean_signal_coarse = a0 * sin(2*pi*f0*t_coarse);
    
    % 从细粒度布朗运动路径 dW_ref 构造 coarse 粒度的 dW
    dW_coarse = zeros(n_steps_coarse, 1);
    for k = 1:n_steps_coarse
        idx_start = (k-1)*m + 1;
        idx_end   = k*m;
        dW_coarse(k) = sum(dW_ref(idx_start:idx_end));
    end
    
    % 构造 coarse 粒度的随机力 noise_seq_coarse，使：
    %   dx_stoch = noise_seq_coarse(k) * h_coarse = sqrt(2D) dW_coarse(k)
    noise_seq_coarse = zeros(n_coarse, 1);
    noise_seq_coarse(1:n_steps_coarse) = sqrt(2 * d_noise) * dW_coarse / h_coarse;
    noise_seq_coarse(end) = 0;
    
    % 使用 RK4Solver2 在较粗采样率下计算轨迹
    x_coarse_hsubsr = RK4Solver2(drift_hsubsr, clean_signal_coarse, noise_seq_coarse, fs_coarse);
    x_coarse_ubsr = RK4Solver2(drift_ubsr, clean_signal_coarse, noise_seq_coarse, fs_coarse);
    x_coarse_plbsr  = RK4Solver2(drift_plbsr,  clean_signal_coarse, noise_seq_coarse, fs_coarse);
    
    % 从参考解中抽取与粗栅格对应的采样点 (保持时间对齐)
    % 参考解的索引: 1, 1+m, 1+2m, ..., 1+n_steps_coarse*m = n_ref
    idx_ref = 1:m:n_ref;
    x_ref_hsubsr_sub = x_ref_hsubsr(idx_ref);
    x_ref_ubsr_sub = x_ref_ubsr(idx_ref);
    x_ref_plbsr_sub  = x_ref_plbsr(idx_ref);
    
    % 计算 RMS 轨迹误差
    rms_err_hsubsr(i_fs) = sqrt(mean((x_coarse_hsubsr - x_ref_hsubsr_sub).^2));
    rms_err_ubsr(i_fs) = sqrt(mean((x_coarse_ubsr - x_ref_ubsr_sub).^2));
    rms_err_plbsr(i_fs)  = sqrt(mean((x_coarse_plbsr  - x_ref_plbsr_sub ).^2));
    
    fprintf('[fs = %4d Hz] RMS误差: HSUBSR = %.4f, UBSR = %.4f, PLBSR = %.4f\n', ...
        fs_coarse, rms_err_hsubsr(i_fs), rms_err_ubsr(i_fs), rms_err_plbsr(i_fs));
end

%% 5. 绘制 RMS 误差 - fs 曲线 ==========================================

SetThesisDefaultStyle();
fig = CreateThesisFigure(8, 6, 2);

subplot(2,1,1);
plot(fs_list, rms_err_hsubsr, 'o--', 'LineWidth', 1.5); hold on;
plot(fs_list, rms_err_ubsr, 's-.', 'LineWidth', 1.5);
plot(fs_list, rms_err_plbsr,  'd:', 'LineWidth', 1.5);
xlabel('$f_s$', 'Interpreter', 'latex');
ylabel('RMS Error');
title('Comparison of Trajectory Errors');
legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'northeast');
grid on;

subplot(2,1,2);
semilogy(fs_list, rms_err_hsubsr + eps, 'o--', 'LineWidth', 1.5); hold on;
semilogy(fs_list, rms_err_ubsr + eps, 's-.', 'LineWidth', 1.5);
semilogy(fs_list, rms_err_plbsr  + eps, 'd:', 'LineWidth', 1.5);
xlabel('$f_s$', 'Interpreter', 'latex');
ylabel('RMS Error (Log Scale)');
title('Comparison of Trajectory Errors (Log Scale)');
legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'northeast');
grid on;