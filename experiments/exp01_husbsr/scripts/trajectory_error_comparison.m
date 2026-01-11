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

%% 1. 全局仿真参数 =====================================================

t_total      = 1000;     % 总仿真时间 (s)
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

%% 2. 构造参考时间栅格与确定性信号 =====================================

h_ref        = 1 / fs_ref;
n_steps_ref  = t_total * fs_ref;          % 参考步数
n_ref        = n_steps_ref + 1;           % 参考采样点数
t_ref        = (0:n_ref-1)' * h_ref;

% 参考信号 (高采样率)
clean_signal_ref = a0 * sin(2*pi*f0*t_ref);

%% 3. 蒙特卡洛仿真 (并行计算) ==========================================

n_mc = 50; % 蒙特卡洛仿真次数
fprintf('Starting Monte Carlo Simulation with %d trials...\n', n_mc);

% 检查采样率整除关系
for i_fs = 1:num_fs
    if mod(fs_ref, fs_list(i_fs)) ~= 0
        error('fs_ref (%d) 必须能被 fs_list 中的每个采样率 (%d) 整除！', fs_ref, fs_list(i_fs));
    end
end

% 预分配空间存储每次实验的误差结果
rms_err_hsubsr_trials = zeros(num_fs, n_mc);
rms_err_ubsr_trials   = zeros(num_fs, n_mc);
rms_err_plbsr_trials  = zeros(num_fs, n_mc);

% 定义 drift 函数 (broadcast 到 workers)
drift_hsubsr = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);
drift_ubsr   = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
drift_plbsr  = @(x) PLBSR_Dynamics(x, u_plbsr, l_plbsr);

% 开启并行池 (如果尚未开启)
if isempty(gcp('nocreate'))
    parpool;
end

% 并行循环
parfor i_mc = 1:n_mc
    % 3.1 生成随机噪声 (dW_ref 对于每次 Trial 是独立的)
    % dW_ref(i) ~ N(0, h_ref)
    dW_ref = sqrt(h_ref) * randn(n_steps_ref, 1);
    
    noise_seq_ref = zeros(n_ref, 1);
    noise_seq_ref(1:n_steps_ref) = sqrt(2 * d_noise) * dW_ref / h_ref;
    noise_seq_ref(end) = 0; % 占位
    
    % 3.2 计算参考轨迹 (Ground Truth)
    x_ref_hsubsr_mc = RK4Solver2(drift_hsubsr, clean_signal_ref, noise_seq_ref, fs_ref);
    x_ref_ubsr_mc   = RK4Solver2(drift_ubsr, clean_signal_ref, noise_seq_ref, fs_ref);
    x_ref_plbsr_mc  = RK4Solver2(drift_plbsr,  clean_signal_ref, noise_seq_ref, fs_ref);
    
    % 临时变量存储单次实验的所有 fs 结果
    tmp_err_hsubsr = zeros(num_fs, 1);
    tmp_err_ubsr   = zeros(num_fs, 1);
    tmp_err_plbsr  = zeros(num_fs, 1);
    
    % 3.3 遍历不同粗采样率
    for i_fs = 1:num_fs
        fs_coarse = fs_list(i_fs);
        h_coarse  = 1 / fs_coarse;
        m         = fs_ref / fs_coarse;
        
        n_steps_coarse = t_total * fs_coarse;
        n_coarse       = n_steps_coarse + 1;
        
        t_coarse            = (0:n_coarse-1)' * h_coarse;
        clean_signal_coarse = a0 * sin(2*pi*f0*t_coarse);
        
        % 降采样噪声 (累加 Wiener 增量)
        dW_coarse = zeros(n_steps_coarse, 1);
        for k = 1:n_steps_coarse
            idx_start = (k-1)*m + 1;
            idx_end   = k*m;
            dW_coarse(k) = sum(dW_ref(idx_start:idx_end));
        end
        
        noise_seq_coarse = zeros(n_coarse, 1);
        noise_seq_coarse(1:n_steps_coarse) = sqrt(2 * d_noise) * dW_coarse / h_coarse;
        noise_seq_coarse(end) = 0;
        
        % 计算粗粒度轨迹
        x_coarse_h = RK4Solver2(drift_hsubsr, clean_signal_coarse, noise_seq_coarse, fs_coarse);
        x_coarse_u = RK4Solver2(drift_ubsr, clean_signal_coarse, noise_seq_coarse, fs_coarse);
        x_coarse_p = RK4Solver2(drift_plbsr,  clean_signal_coarse, noise_seq_coarse, fs_coarse);
        
        % 降采样参考轨迹以对齐
        idx_ref = 1:m:n_ref;
        x_ref_h_sub = x_ref_hsubsr_mc(idx_ref);
        x_ref_u_sub = x_ref_ubsr_mc(idx_ref);
        x_ref_p_sub = x_ref_plbsr_mc(idx_ref);
        
        % 计算误差
        tmp_err_hsubsr(i_fs) = sqrt(mean((x_coarse_h - x_ref_h_sub).^2));
        tmp_err_ubsr(i_fs)   = sqrt(mean((x_coarse_u - x_ref_u_sub).^2));
        tmp_err_plbsr(i_fs)  = sqrt(mean((x_coarse_p - x_ref_p_sub).^2));
    end
    
    % 存储该次实验结果
    rms_err_hsubsr_trials(:, i_mc) = tmp_err_hsubsr;
    rms_err_ubsr_trials(:, i_mc)   = tmp_err_ubsr;
    rms_err_plbsr_trials(:, i_mc)  = tmp_err_plbsr;
    
    % parfor 中不建议频繁使用 fprintf，但偶尔使用用于调试可以
    % fprintf('Finished trial %d\n', i_mc);
end

%% 4. 统计结果 =========================================================

rms_err_hsubsr = mean(rms_err_hsubsr_trials, 2);
rms_err_ubsr   = mean(rms_err_ubsr_trials, 2);
rms_err_plbsr  = mean(rms_err_plbsr_trials, 2);

% 打印结果 (Mean)
fprintf('Completed %d Monte Carlo trials.\n', n_mc);
for i_fs = 1:num_fs
    fprintf('[fs = %4d Hz] Mean RMS Error: HSUBSR = %.4f, UBSR = %.4f, PLBSR = %.4f\n', ...
        fs_list(i_fs), rms_err_hsubsr(i_fs), rms_err_ubsr(i_fs), rms_err_plbsr(i_fs));
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

subplot(2,1,2);
semilogy(fs_list, rms_err_hsubsr + eps, 'o--', 'LineWidth', 1.5); hold on;
semilogy(fs_list, rms_err_ubsr + eps, 's-.', 'LineWidth', 1.5);
semilogy(fs_list, rms_err_plbsr  + eps, 'd:', 'LineWidth', 1.5);
xlabel('$f_s$', 'Interpreter', 'latex');
ylabel('RMS Error (Log Scale)');
title('Comparison of Trajectory Errors (Log Scale)');
legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'northeast');

%% 6. 保存结果 =========================================================
results.fs_list = fs_list;
results.rms_err_hsubsr = rms_err_hsubsr;
results.rms_err_ubsr   = rms_err_ubsr;
results.rms_err_plbsr  = rms_err_plbsr;