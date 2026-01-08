% =========================================================================
% Description: 随机采样 HSUBSR 势函数参数，在每个样本上计算：
%              - 盲适应度指标 RSCM / J_total
%              - 盲共振度 eta / η
%              - 基于真实 f0 的输出 SNRo
%              分析三者之间的相关性，用于验证盲指标和盲共振度的合理性。
%
% Author: LiuShuang
% Created: 2026-01-06
% Last Modified: 2026-01-06
% =========================================================================

clc; clear; close all;

%% 1. 公共仿真参数设置 ==============================================
fs         = 5;                 % 采样频率 Hz
T          = 2000;              % 信号时长 s
n_samples  = fs * T;
t_vec      = (0:n_samples-1)' / fs;

A0         = 0.05;              % 弱信号幅值
f0    = 0.01;                   % 真实信号频率，仅用于 SNR 验证
clean_sig  = A0 * sin(2*pi*f0*t_vec);

noise_intensity = 0.1;          % 噪声强度 D
steady_ratio    = 0.1;          % 丢弃瞬态比例
h_step          = 1/fs;
noise_seq = sqrt(2*noise_intensity/h_step) * randn(n_samples, 1);

% HSUBSR 参数范围
theta_min = [0.5,  0.5, 0.5, 0.5];
theta_max = [5.0,  5.0, 5.0, 5.0];

opts.kappa_eta = 5.0;     % 控制 eta 的 sigmoid 陡峭性
opts.w_ami     = 0.1;     % 适应度权重
opts.w_mpe     = 0.9;

%% 2. 随机采样参数并计算 (J, eta, SNR) ================================
num_samples_param = 100;   % 采样个数，可根据时间调整

theta_dim   = numel(theta_min);

J_total_list   = zeros(num_samples_param, 1);
eta_list       = zeros(num_samples_param, 1);
snr_list       = zeros(num_samples_param, 1);

J_occ_list     = zeros(num_samples_param, 1);
J_ami_list     = zeros(num_samples_param, 1);
J_entropy_list = zeros(num_samples_param, 1);

for i = 1:num_samples_param
    theta_rand = theta_min + rand(1, theta_dim) .* (theta_max - theta_min);
    
    % ---- 计算盲适应度 & 盲共振度 ----
    [J_fit, eta, dbg] = RSCMEvaluator(theta_rand, clean_sig, noise_seq, fs, opts);
    
    J_total_list(i)   = J_fit;
    eta_list(i)       = eta;
    
    J_occ_list(i)     = dbg.J_occ;
    J_ami_list(i)     = dbg.J_ami;
    J_entropy_list(i) = dbg.J_mpe;
    
    % ---- 计算对应的 SNR (使用真实 f0，仅用于分析) ----
    drift = @(x) HSUBSR_Dynamics(x, theta_rand(1), theta_rand(2), ...
        theta_rand(3), theta_rand(4));
    x_out = RK4Solver2(drift, clean_sig, noise_seq, fs);
    x_out = x_out(round(steady_ratio*n_samples):end);
    
    snr_list(i) = SNRo2(x_out, fs, f0);
end

%% 4. 计算相关系数并打印 ================================================
corr_eta_J      = corr(eta_list, J_total_list);
corr_eta_snr    = corr(eta_list, snr_list);
corr_J_snr      = corr(J_total_list, snr_list);

fprintf('corr(eta, J_total)   = %.3f\n', corr_eta_J);
fprintf('corr(eta, SNR)       = %.3f\n', corr_eta_snr);
fprintf('corr(J_total, SNR)   = %.3f\n', corr_J_snr);

%% 5. 可视化 eta 分布 & eta–J/SNR 散点图 ================================
figure('Color','w');
subplot(2,2,1);
histogram(eta_list, 10);
xlabel('\eta'); ylabel('计数');
title('\eta 分布');
grid on;

subplot(2,2,2);
scatter(eta_list, J_total_list, 'filled');
xlabel('\eta'); ylabel('J_{total}');
title(sprintf('\n corr(eta, J_{total}) = %.3f', corr_eta_J));
grid on;

subplot(2,2,3);
scatter(eta_list, snr_list, 'filled');
xlabel('\eta'); ylabel('SNR (dB)');
title(sprintf('corr(eta, SNR) = %.3f', corr_eta_snr));
grid on;

subplot(2,2,4);
scatter(J_total_list, snr_list, 'filled');
xlabel('J_{total}'); ylabel('SNR (dB)');
title(sprintf('corr(J_{total}, SNR) = %.3f', corr_J_snr));
grid on;