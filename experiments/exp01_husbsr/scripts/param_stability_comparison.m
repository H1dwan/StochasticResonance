% =========================================================================
% Description: 参数扰动鲁棒性对比
%              (HSUBSR vs UBSR\PLBSR)
%
% Author: LiuShuang
% Created: 2025-12-04
% Last Modified: 2025-12-22
%
% Usage:
%   当前存在一定问题，需要重新调试
% =========================================================================

clc; clear; close all;

%% 1. 全局仿真参数设置 ===============================================

% --------- 信号与噪声参数 -----------------------------------------
T_total         = 1000;       % 仿真总时长 (s)，越大统计量越稳定，可视内存情况调整
f0              = 0.01;       % 弱周期信号频率 (Hz)，满足绝热近似
A0              = 0.1;        % 弱周期信号幅值 (远小于势垒高度)
D0              = 1.0;        % 用于仿真4的固定噪声强度（接近SR最优噪声，可根据实验微调）
TRANSIENT_RATIO = 0.1;        % 前 10% 作为瞬态丢弃

% --------- 步长/采样率列表 (仿真4A) -------------------------------
fs_list         = 5:5:50;   % 采样频率列表 (Hz)，fs越大步长越小
num_fs          = numel(fs_list);

% -------- 固定的势函数结构参数 ----------------------------------------
xm = 1;
dU = 0.25;

% -------- HSUBSR 势函数参数 ------------------------------------------
shape_factor = 8;
[a_hsubsr0, b_hsubsr0, k1_hsubsr0, k2_hsubsr0] = CalibrateHSUBSR(xm, dU, shape_factor);

% -------- UBSR 势函数参数 --------------------------------------------
a_ubsr0     = 1.0;
b_ubsr0     = 1.0;

% -------- PLBSR 势函数参数 -------------------------------------------
u_plbsr0    = 0.25;
l_plbsr0    = 1;

% --------- Monte Carlo 重复次数 ------------------------------------
n_repeat_4A = 10;   % 仿真4A中每个 fs 的重复次数
n_repeat_4B = 12;   % 仿真4B中每组参数扰动的样本数

%% 2. 仿真4A：步长 / 采样率敏感性分析 ================================

% % 预分配结果数组
% snr_hsubsr_fs  = zeros(num_fs, 1);
% snr_ubsr_fs  = zeros(num_fs, 1);
% snr_plbsr_fs  = zeros(num_fs, 1);

% for i_fs = 1:num_fs
%     fs = fs_list(i_fs);
%     h  = 1 / fs;                      % 积分步长

%     N  = round(T_total * fs);        % 采样点数
%     t  = (0:N-1)' / fs;              % 时间向量
%     clean_signal = A0 * sin(2*pi*f0*t);  % 弱周期信号

%     idx0 = round(TRANSIENT_RATIO * N);     % 瞬态截断位置

%     % 用于累积 Monte Carlo 的 SNR
%     snr_hsubsr_rep = zeros(n_repeat_4A, 1);
%     snr_ubsr_rep = zeros(n_repeat_4A, 1);
%     snr_plbsr_rep  = zeros(n_repeat_4A, 1);

%     for k = 1:n_repeat_4A
%         % 生成噪声：xi(t)，满足 int xi dt ≈ noise_seq * h
%         % 对应随机微分方程 dx = ... + sqrt(2D) dW
%         % 令 noise_seq = sqrt(2D/h)*randn → dx_stoch ≈ noise_seq*h = sqrt(2D*h)*randn
%         noise_seq = sqrt(2 * D0 * fs) * randn(N, 1);

%         % ----------------- HESUBSR 模型 -------------------
%         drift_hsubsr = @(x) HSUBSR_Dynamics(x, a_hsubsr0, b_hsubsr0, k1_hsubsr0, k2_hsubsr0);
%         x_hsubsr = RK4Solver2(drift_hsubsr, clean_signal, noise_seq, fs);
%         x_hsubsr_stead = x_hsubsr(idx0:end);

%         % 计算输出 SNR
%         snr_hsubsr_rep(k) = SNRo(x_hsubsr_stead, fs, f0);

%         % ----------------- UBSR 模型 ---------------------
%         drift_ubsr = @(x) UBSR_Dynamics(x, a_ubsr0, b_ubsr0);
%         x_ubsr = RK4Solver2(drift_ubsr, clean_signal, noise_seq, fs);
%         x_ubsr_stead = x_ubsr(idx0:end);

%         snr_ubsr_rep(k) = SNRo(x_ubsr_stead, fs, f0);

%         % ----------------- PLBSR 模型 ---------------------
%         drift_plbsr = @(x) PLBSR_Dynamics(x, u_plbsr0, l_plbsr0);
%         x_plbsr = RK4Solver2(drift_plbsr, clean_signal, noise_seq, fs);
%         x_plbsr_stead = x_plbsr(idx0:end);

%         snr_plbsr_rep(k) = SNRo(x_plbsr_stead, fs, f0);

%     end

%     % 对 Monte Carlo 重复求平均
%     snr_hsubsr_fs(i_fs) = mean(snr_hsubsr_rep);
%     snr_ubsr_fs(i_fs) = mean(snr_ubsr_rep);
%     snr_plbsr_fs(i_fs) = mean(snr_plbsr_rep);

%     fprintf('[4A] fs = %4d Hz  HESUBSR-SNR = %.3f dB, UBSR-SNR = %.3f dB\n', ...
%         fs, snr_hsubsr_fs(i_fs), snr_ubsr_fs(i_fs));
% end

% % --- 以最高采样率 fs_max 结果作为“近似真值”，计算相对误差 ----------
% snr_hsubsr_ref = snr_hsubsr_fs(end);
% snr_ubsr_ref = snr_ubsr_fs(end);
% snr_plbsr_ref = snr_plbsr_fs(end);

% relerr_snr_hsubsr = abs(snr_hsubsr_fs - snr_hsubsr_ref) ./ abs(snr_hsubsr_ref);
% relerr_snr_ubsr = abs(snr_ubsr_fs - snr_ubsr_ref) ./ abs(snr_ubsr_ref);
% relerr_snr_plbsr = abs(snr_plbsr_fs - snr_plbsr_ref) ./ abs(snr_plbsr_ref);
%% 2.1 绘制 SNR vs fs 曲线 ==========================================

% figure('Position', [100, 100, 1200, 450], 'Color', 'white');
% plot(fs_list, snr_hsubsr_fs, 'o-', 'LineWidth', 1.5); hold on;
% plot(fs_list, snr_ubsr_fs, 's--', 'LineWidth', 1.5);
% plot(fs_list, snr_plbsr_fs, 'x:', 'LineWidth', 1.5);
% xlabel('采样频率 f_s (Hz)');
% ylabel('输出 SNR (dB)');
% title('不同采样率下输出 SNR 对比');
% legend({'HSUBSR','UBSR', 'PLBSR'}, 'Location', 'best');
% grid on;

%% 2.2 绘制相对误差 vs fs 曲线 ======================================

% figure('Position', [100, 600, 1200, 450], 'Color', 'white');

% semilogy(fs_list, relerr_snr_hsubsr + eps, 'o-', 'LineWidth', 1.5); hold on;
% semilogy(fs_list, relerr_snr_ubsr + eps, 's--', 'LineWidth', 1.5);
% semilogy(fs_list, relerr_snr_plbsr + eps, 'x:', 'LineWidth', 1.5);
% xlabel('采样频率 f_s (Hz)');
% ylabel('SNR 相对误差');
% title('SNR 对采样率的敏感性');
% legend({'HSUBSR','UBSR', 'PLBSR'}, 'Location', 'best');
% grid on;

%% 3. 仿真4B：参数扰动鲁棒性分析 =====================================

fs_robust = 200;              % 参数扰动仿真时固定采样频率
h_robust  = 1 / fs_robust;
N_robust  = round(T_total * fs_robust);
t_robust  = (0:N_robust-1)' / fs_robust;
clean_signal_robust = A0 * sin(2*pi*f0 * t_robust);
idx0_robust = round(TRANSIENT_RATIO * N_robust);

% 参数扰动幅度（相对误差）
delta_max = 0.05;  % ±5%

snr_hsubsr_samples = zeros(n_repeat_4B, 1);
snr_ubsr_samples = zeros(n_repeat_4B, 1);
snr_plbsr_samples = zeros(n_repeat_4B, 1);

parfor i_sample = 1:n_repeat_4B
    % -------- 生成一条公共噪声轨迹，用于对比两种势 ---------------
    noise_seq_robust = sqrt(2 * D0 / h_robust) * randn(N_robust, 1);
    
    % -------- HSUBSR 参数扰动 -----------------------------------
    a_hsubsr = a_hsubsr0 * (1 + delta_max * (2*rand - 1));
    b_hsubsr = b_hsubsr0 * (1 + delta_max * (2*rand - 1));
    k1_hsubsr = k1_hsubsr0 * (1 + delta_max * (2*rand - 1));
    k2_hsubsr = k2_hsubsr0 * (1 + delta_max * (2*rand - 1));
    
    drift_hsubsr_pert = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);
    x_hsubsr_pert = RK4Solver2(drift_hsubsr_pert, clean_signal_robust, noise_seq_robust, fs_robust);
    x_hsubsr_stead = x_hsubsr_pert(idx0_robust:end);
    snr_hsubsr_samples(i_sample) = SNRo2(x_hsubsr_stead, fs_robust, f0);
    
    % -------- UBSR 参数扰动 -------------------------------------
    a_ubsr = a_ubsr0 * (1 + delta_max * (2*rand - 1));
    b_ubsr = b_ubsr0 * (1 + delta_max * (2*rand - 1));
    
    drift_ubsr_pert = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
    x_ubsr_pert = RK4Solver2(drift_ubsr_pert, clean_signal_robust, noise_seq_robust, fs_robust);
    x_ubsr_stead = x_ubsr_pert(idx0_robust:end);
    snr_ubsr_samples(i_sample) = SNRo2(x_ubsr_stead, fs_robust, f0);
    
    % -------- PLBSR 模型参数扰动 ----------------------------------
    u_plbsr = u_plbsr0 * (1 + delta_max * (2*rand - 1));
    l_plbsr = l_plbsr0 * (1 + delta_max * (2*rand - 1));
    
    drift_plbsr_pert = @(x) PLBSR_Dynamics(x, u_plbsr, l_plbsr);
    x_plbsr_pert = RK4Solver2(drift_plbsr_pert, clean_signal_robust, noise_seq_robust, fs_robust);
    x_plbsr_stead = x_plbsr_pert(idx0_robust:end);
    snr_plbsr_samples(i_sample) = SNRo2(x_plbsr_stead, fs_robust, f0);
    
end

% 计算统计量：均值和标准差
snr_hsubsr_mean = mean(snr_hsubsr_samples);
snr_hsubsr_std  = std(snr_hsubsr_samples);
snr_ubsr_mean = mean(snr_ubsr_samples);
snr_ubsr_std  = std(snr_ubsr_samples);
snr_plbsr_mean = mean(snr_plbsr_samples);
snr_plbsr_std  = std(snr_plbsr_samples);

fprintf('\n[4B] HSUBSR 参数扰动: SNR 均值 = %.3f dB, 标准差 = %.3f dB\n', ...
    snr_hsubsr_mean, snr_hsubsr_std);
fprintf('[4B] UBSR 参数扰动:   SNR 均值 = %.3f dB, 标准差 = %.3f dB\n', ...
    snr_ubsr_mean, snr_ubsr_std);
fprintf('[4B] PLBSR 模型参数扰动: SNR 均值 = %.3f dB, 标准差 = %.3f dB\n', ...
    snr_plbsr_mean, snr_plbsr_std);

%% 3.1 绘制参数扰动下 SNR 分布 =======================================

figure('Position', [400, 200, 900, 500], 'Color', 'white');

subplot(1,2,1);
boxplot([snr_hsubsr_samples, snr_ubsr_samples, snr_plbsr_samples], ...
    'Labels', {'HSUBSR','UBSR', 'PLBSR'});
ylabel('输出 SNR (dB)');
title(sprintf('参数扰动 (±%.1f%%) 下 SNR 分布', delta_max*100));
grid on;

subplot(1,2,2);
histogram(snr_hsubsr_samples, 10); hold on;
histogram(snr_ubsr_samples, 10);
histogram(snr_plbsr_samples, 10);
xlabel('输出 SNR (dB)');
ylabel('出现次数');
legend({'HSUBSR','UBSR', 'PLBSR'}, 'Location', 'best');
title('参数扰动下 SNR 直方图对比');
grid on;