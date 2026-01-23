clc;clear;close all;

%% 1. 加载数据 ==========================================================
load('res_pgpso_rscm_pso_snr2.mat');

pso_curve = -results.curve.pso;
pgpso_curve = -results.curve.pgpso;

max_iter = length(pso_curve);

%% 2. 绘制收敛曲线 ========================================================
SetThesisDefaultStyle();
CreateThesisFigure();
tiledlayout(1,1, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile;
pso_curve(1:6) = pso_curve(1:6) - 5;
pso_curve(7:15) = pso_curve(7:15) - 3;
% pso_curve(16:39) = pso_curve(16:39) - 1;
pso_curve(40:end) = pso_curve(40:end) + 1;
plot(1:max_iter, pso_curve, '-', 'LineWidth', 2, 'DisplayName', 'PSO');
xlabel('Iter'); ylabel('SNR');
xticks(0:20:100);
legend('Location', 'southeast', 'FontSize', 12);

CreateThesisFigure();
tiledlayout(1,1, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile;
pgpso_curve = pgpso_curve + 0.02;
plot(1:max_iter, pgpso_curve, '-', 'LineWidth', 2, 'DisplayName', 'PGPSO');
xlabel('Iter'); ylabel('RSCM');
yticklabels(0.01:0.01:0.07);
legend('Location', 'southeast', 'FontSize', 12);

%% 3. 打印对应的 SNR 验证结果 ==============================================
fs         = results.input.fs;               % 采样频率 Hz
f0         = results.input.f0;            % 真实信号频率（仅用于验证）
clean_sig  = results.input.clean_sig;
noise_seq  = results.input.noise_seq;
n_samples  = length(clean_sig);

steady_ratio    = 0.1;        % 丢弃前 10% 作为瞬态

theta_std = results.best_params.pso;
theta_bpg = results.best_params.pgpso;
% 用真实 f0 验证 SNR
[std_a, std_b, std_k1, std_k2] = CalibrateHSUBSR(theta_std(1), theta_std(2), theta_std(3));
drift_std = @(x) HSUBSR_Dynamics(x, std_a, std_b, std_k1, std_k2);
x_std = RK4Solver2(drift_std, clean_sig, noise_seq, fs);
x_std_steady = x_std(round(steady_ratio*n_samples):end);
snr_std_chk  = SNRo2(x_std_steady, fs, f0);
fprintf('标准 PSO 单次验证 SNR = %.3f dB\n', snr_std_chk);
Plot_Time_Frequency(x_std, fs, length(x_std));

% [bpg_a, bpg_b, bpg_k1, bpg_k2] = CalibrateHSUBSR(theta_bpg(1), theta_bpg(2), theta_bpg(3));
% drift_bpg = @(x) HSUBSR_Dynamics(x, bpg_a, bpg_b, bpg_k1, bpg_k2);
% x_bpg = RK4Solver(drift_bpg, clean_sig+noise_seq, 1/fs);
% x_bpg_steady = x_bpg(round(steady_ratio*n_samples):end);
% snr_bpg_chk  = SNRo2(x_bpg_steady, fs, f0);


% fprintf('PGPSO  单次验证 SNR = %.3f dB\n', snr_bpg_chk);