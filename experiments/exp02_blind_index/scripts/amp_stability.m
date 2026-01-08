% =========================================================================
% Description: 验证 RSCM 指标的幅值稳定性
%
% Author: LiuShuang
% Created: 2025-12-25
% Last Modified: 2025-12-25
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 仿真全局参数设置 ===============================================
fs          = 5;               % 采样频率 (Hz)
t_total     = 2000;             % 总仿真时间长度 (秒)
num_samples = fs * t_total;     % 样本点数
time_axis   = (0:num_samples-1)' / fs;

% 输入信号参数（弱周期信号）
f0      = 0.01;               % 目标信号频率 (Hz)，满足 f0 << fs
% A0      = 0.05;               % 信号幅值（可根据需要调整）
% amp_list = [0.05, 0.10, 0.15];
amp_list = [0.15];  % 待测试频率列表
num_amp = length(amp_list);

% 噪声强度 D 扫描范围
D_list   = 0.05:0.01:0.45;
num_D    = length(D_list);

% 每个 D 重复次数（用于统计平均）
n_repeat = 12;

% 结果存储
results = struct('A0', cell(num_amp, 1), ...
    'D_axis', cell(num_amp, 1), ...
    'snr_mean', cell(num_amp, 1), ...
    'ami_mean', cell(num_amp, 1), ...
    'mpe_mean', cell(num_amp, 1), ...
    'rscm_mean', cell(num_amp, 1), ...
    'rho_rscm', cell(num_amp, 1), ...
    'rho_ami', cell(num_amp, 1), ...
    'rho_mpe', cell(num_amp, 1), ...
    'idx_peak_snr', cell(num_amp, 1), ...
    'idx_peak_rscm', cell(num_amp, 1), ...
    'idx_peak_ami', cell(num_amp, 1), ...
    'idx_peak_mpe', cell(num_amp, 1), ...
    'D_peak_snr', cell(num_amp, 1), ...
    'D_peak_rscm', cell(num_amp, 1), ...
    'D_peak_ami', cell(num_amp, 1), ...
    'D_peak_mpe', cell(num_amp, 1), ...
    'delta_peak_rscm', cell(num_amp, 1), ...
    'delta_peak_ami', cell(num_amp, 1), ...
    'delta_peak_mpe', cell(num_amp, 1));
%% 2. 定义双稳系统漂移函数 ==========================================
% drift_func = @(x) CBSR_Dynamics(x, 1, 1);
% drift_func = @(x) UBSR_Dynamics(x, 1, 1);
[a, b, k1, k2] = CalibrateHSUBSR(1, 0.25, 50);
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);

%% 3. 频率稳定性验证 ==========================
for idx_a = 1 : num_amp
    A0 = amp_list(idx_a);
    fprintf('正在仿真幅值 A0 = %.3f (%d/%d)...\n', A0, idx_a, num_amp);
    clean_signal = A0 * sin(2*pi*f0*time_axis);
    
    snr_mean  = zeros(num_D, 1);
    ami_mean  = zeros(num_D, 1);
    mpe_mean  = zeros(num_D, 1);
    rscm_mean = zeros(num_D, 1);
    
    for idx_D = 1:num_D
        D = D_list(idx_D);
        
        snr_rep  = zeros(n_repeat, 1);
        ami_rep  = zeros(n_repeat, 1);
        mpe_rep  = zeros(n_repeat, 1);
        rscm_rep = zeros(n_repeat, 1);
        
        parfor k = 1:n_repeat
            noise_seq = sqrt(2*D*fs) * randn(num_samples, 1);
            x = RK4Solver2(drift_func, clean_signal, noise_seq, fs);
            
            steady_start = round(0.1 * num_samples);
            x_stead = x(steady_start+1:end);
            
            snr_rep(k) = SNRo2(x_stead, fs, f0);
            [ami_rep(k), ~] = AMI2(x_stead, fs);
            
            [scales, ~] = SelectMpeScales(x_stead, fs, 3);
            [~, mpnorm, ~, ~] = MultiScalePermEn(x_stead, scales);
            mpe_val = std(mpnorm(~isnan(mpnorm)));
            mpe_rep(k) = mpe_val;
            
            rscm_rep(k) = ami_rep(k)^0.5 * (1 - mpe_val)^0.5;
        end
        
        snr_mean(idx_D)  = mean(snr_rep);
        ami_mean(idx_D)  = mean(ami_rep);
        mpe_mean(idx_D)  = mean(mpe_rep);
        rscm_mean(idx_D) = mean(rscm_rep);
    end
    
    [~, idx_peak_snr]  = max(snr_mean);
    [~, idx_peak_rscm] = max(rscm_mean);
    [~, idx_peak_ami]  = max(ami_mean);
    [~, idx_peak_mpe]  = max(1 - mpe_mean);
    
    D_peak_snr   = D_list(idx_peak_snr);
    D_peak_rscm  = D_list(idx_peak_rscm);
    D_peak_ami   = D_list(idx_peak_ami);
    D_peak_mpe   = D_list(idx_peak_mpe);
    
    R_rscm = corrcoef(snr_mean, rscm_mean);
    R_ami  = corrcoef(snr_mean, ami_mean);
    R_mpe  = corrcoef(snr_mean, 1 - mpe_mean);
    
    results(idx_a).A0 = A0;
    results(idx_a).D_axis = D_list;
    results(idx_a).snr_mean = snr_mean;
    results(idx_a).ami_mean = ami_mean;
    results(idx_a).mpe_mean = mpe_mean;
    results(idx_a).rscm_mean = rscm_mean;
    results(idx_a).rho_rscm = R_rscm(1, 2);
    results(idx_a).rho_ami = R_ami(1, 2);
    results(idx_a).rho_mpe = R_mpe(1, 2);
    results(idx_a).idx_peak_snr = idx_peak_snr;
    results(idx_a).idx_peak_rscm = idx_peak_rscm;
    results(idx_a).idx_peak_ami = idx_peak_ami;
    results(idx_a).idx_peak_mpe = idx_peak_mpe;
    results(idx_a).D_peak_snr = D_peak_snr;
    results(idx_a).D_peak_rscm = D_peak_rscm;
    results(idx_a).D_peak_ami = D_peak_ami;
    results(idx_a).D_peak_mpe = D_peak_mpe;
    results(idx_a).delta_peak_rscm = abs(D_peak_rscm - D_peak_snr);
    results(idx_a).delta_peak_ami = abs(D_peak_ami - D_peak_snr);
    results(idx_a).delta_peak_mpe = abs(D_peak_mpe - D_peak_snr);
end

%% 4. 结果打印
fprintf('\n================ 频率稳定性结果 ================\n');
for idx_a = 1:num_amp
    fprintf('A0=%.3f | D_peak(SNR)=%.4f | D_peak(RSCM)=%.4f | ΔD=%.4f | corr=%.3f\n', ...
        results(idx_a).A0, results(idx_a).D_peak_snr, results(idx_a).D_peak_rscm, ...
        results(idx_a).delta_peak_rscm, results(idx_a).rho_rscm);
    fprintf('          | D_peak(AMI)=%.4f | ΔD=%.4f | corr=%.3f\n', ...
        results(idx_a).D_peak_ami, results(idx_a).delta_peak_ami, results(idx_a).rho_ami);
    fprintf('          | D_peak(MPE)=%.4f | ΔD=%.4f | corr=%.3f\n', ...
        results(idx_a).D_peak_mpe, results(idx_a).delta_peak_mpe, results(idx_a).rho_mpe);
end
fprintf('=================================================\n');

%% 5. 曲线绘制（按频率分面）
SetThesisDefaultStyle();
fig1 = CreateThesisFigure();
for idx_a = 1:num_amp
    subplot(2, 2, idx_a);
    D_axis = results(idx_a).D_axis;
    
    plot(D_axis, smooth(results(idx_a).ami_mean, 3), 'o-', 'DisplayName', 'SNR'); hold on;
    % plot(D_axis, smooth(results(idx_a).ami_mean, 1), 'd-.', 'DisplayName', 'AMI');
    % plot(D_axis, smooth(1 - results(idx_a).mpe_mean, 1), 's--', 'DisplayName', '1-MPE');
    % plot(D_axis, smooth(results(idx_a).rscm_mean, 1), '^-', 'DisplayName', 'RSCM');
    
    % plot(results(idx_a).D_peak_snr, results(idx_a).snr_mean(results(idx_a).idx_peak_snr), 'ko', 'MarkerFaceColor', 'k');
    % plot(results(idx_a).D_peak_rscm, results(idx_a).rscm_mean(results(idx_a).idx_peak_rscm), 'kp', 'MarkerFaceColor', 'y');
    
    % title(sprintf('A0 = %.3f, corr = %.3f', results(idx_a).A0, results(idx_a).rho_rscm));
    % xlabel('D');
    ylabel('Metric value');
    legend('Location', 'best');
end

%% 6. 稳定性汇总随频率变化
% fig2 = CreateThesisFigure();
% subplot(1, 2, 1);
% plot(amp_list, [results.D_peak_snr], 'o-', 'DisplayName', 'D_{peak} SNR'); hold on;
% plot(amp_list, [results.D_peak_rscm], 'd-.', 'DisplayName', 'D_{peak} RSCM');
% plot(amp_list, [results.D_peak_ami], 's--', 'DisplayName', 'D_{peak} AMI');
% plot(amp_list, [results.D_peak_mpe], 'x:', 'DisplayName', 'D_{peak} 1-MPE');
% xlabel('A0');
% ylabel('D_{peak}');
% legend('Location', 'best');

% subplot(1, 2, 2);
% yyaxis left;
% plot(amp_list, [results.delta_peak_rscm], 'o-', 'DisplayName', '|ΔD| RSCM'); hold on;
% plot(amp_list, [results.delta_peak_ami], 'd-.', 'DisplayName', '|ΔD| AMI');
% plot(amp_list, [results.delta_peak_mpe], 's--', 'DisplayName', '|ΔD| 1-MPE');
% ylabel('|ΔD| relative to SNR');
% yyaxis right;
% plot(amp_list, [results.rho_rscm], '^-', 'DisplayName', 'corr(SNR,RSCM)');
% plot(amp_list, [results.rho_ami], 'x:', 'DisplayName', 'corr(SNR,AMI)');
% plot(amp_list, [results.rho_mpe], 'p-.', 'DisplayName', 'corr(SNR,1-MPE)');
% ylabel('Correlation');
% xlabel('A0');
% legend('Location', 'best');