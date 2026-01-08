% =========================================================================
% Description: 在同一 SR 条件下，验证 PGPSO-RSCM 与 PSO-SNR 的一致性
%              - PSO-SNR：标准 PSO 算法优化 SNR 适应度
%              - PGPSO-RSCM：基于物理引导的 PSO 优化 RSCM 适应度
%
% Author: LiuShuang
% Created: 2026-01-06
% Last Modified: 2026-01-06
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; close all;
% rng(28);

%% 1. 公共仿真参数 =========================================================
fs         = 5;               % 采样频率 Hz
T          = 2000;            % 信号时长 s
n_samples  = T * fs;          % 信号采样点数
t_vec      = (0:n_samples-1)' / fs;

A0         = 0.1;             % 弱信号幅值
f0         = 0.01;            % 真实信号频率（仅用于验证）
% clean_sig  = A0 * sin(2*pi*f0*t_vec);

steady_ratio    = 0.1;        % 丢弃前 10% 作为瞬态A
noise_intensity = 0.12;       % 基准噪声强度 D
% noise_seq = sqrt(2*noise_intensity*fs) * randn(n_samples,1);

% HSUBSR 势结构参数搜索范围，使用 "物理映射搜索策略" (xm, dU, shape_factor)
input_rms = std(clean_sig+noise_seq);
theta_min = [0.5*input_rms, 0.15, 1.1];  % [xm, dU, shape_factor]
theta_max = [2.0*input_rms, 0.5, 50.0];

%% 2. PGPSO / 标准 PSO 公共参数 =============================================
bpg_opts.num_particles = 20;
bpg_opts.max_iter      = 50;
bpg_opts.goal          = 'min';          % 优化目标：'min' 或 'max'

% PSO 参数范围
bpg_opts.w_min   = 0.3;
bpg_opts.w_max   = 0.9;
bpg_opts.c1_min  = 1.0;
bpg_opts.c1_max  = 2.5;
bpg_opts.c2_min  = 1.0;
bpg_opts.c2_max  = 2.5;

% 盲共振度相关参数
bpg_opts.kappa_eta = 5.0;         % 控制 Blind eta 的 sigmoid 陡峭性
bpg_opts.w_ami = 1;
bpg_opts.w_mpe = 1;

bpg_opts.display      = true;     % 实验 1 开日志，实验 2 可改为 false

%% 3. 单次收敛曲线对比（盲适应度） ============================================
fprintf('==== 单次收敛曲线对比 ====\n');

% --- 标准 PSO ---
evaluator_pso = @(x) SNREvaluator(x, clean_sig, noise_seq, fs, f0);
[fit_std, theta_std, curve] = PSO(bpg_opts.num_particles, bpg_opts.max_iter, theta_min, theta_max, ...
    3, evaluator_pso);

% --- PGPSO ---
evaluator_pgpso = @(x) RSCMEvaluator(x, clean_sig, noise_seq, fs, bpg_opts);
[fit_bpg, theta_bpg, hist_bpg] = PGPSO( ...
    theta_min, theta_max, bpg_opts, evaluator_pgpso);

fprintf('PSO 单次最优适应度 SNR = %.4f\n', fit_std);
fprintf('PGPSO  单次最优盲适应度 J = %.4f\n', fit_bpg);

% 用真实 f0 验证 SNR
[std_a, std_b, std_k1, std_k2] = CalibrateHSUBSR(theta_std(1), theta_std(2), theta_std(3));
drift_std = @(x) HSUBSR_Dynamics(x, std_a, std_b, std_k1, std_k2);
x_std = RK4Solver(drift_std, clean_sig+noise_seq, 1/fs);
x_std_steady = x_std(round(steady_ratio*n_samples):end);
snr_std_chk  = SNRo2(x_std_steady, fs, f0);

[bpg_a, bpg_b, bpg_k1, bpg_k2] = CalibrateHSUBSR(theta_bpg(1), theta_bpg(2), theta_bpg(3));
drift_bpg = @(x) HSUBSR_Dynamics(x, bpg_a, bpg_b, bpg_k1, bpg_k2);
x_bpg = RK4Solver(drift_bpg, clean_sig+noise_seq, 1/fs);
x_bpg_steady = x_bpg(round(steady_ratio*n_samples):end);
snr_bpg_chk  = SNRo2(x_bpg_steady, fs, f0);

fprintf('标准 PSO 单次验证 SNR = %.3f dB\n', snr_std_chk);
fprintf('PGPSO  单次验证 SNR = %.3f dB\n', snr_bpg_chk);

%% 4. 绘制收敛曲线 ==========================================================
figure('Color','w'); hold on; box on;
plot(1:bpg_opts.max_iter, -curve, 'LineWidth', 1.5);
xlabel('Iter'); ylabel('SNR');
title('PSO Convergence Curve');

figure('Color','w'); hold on; box on;
plot(1:bpg_opts.max_iter, -hist_bpg.gbest_fit, 'LineWidth', 1.5, 'LineStyle','--');
xlabel('Iter'); ylabel('RSCM');
title('PGPSO Convergence Curve');

%% 5. 保存结果 ==============================================================
results.curve.pso   = curve;
results.curve.pgpso = hist_bpg.gbest_fit;
results.input.fs = fs;
results.input.T = T;
results.input.f0 = f0;
results.input.A0 = A0;
results.input.clean_sig = clean_sig;
results.input.noise_seq = noise_seq;
results.best_params.pso = theta_std;
results.best_params.pgpso = theta_bpg;