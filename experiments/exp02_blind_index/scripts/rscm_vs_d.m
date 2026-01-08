% =========================================================================
% Description: Refined SCM 指标和 SNR 对噪声 D 的依赖关系
%
% Author: LiuShuang
% Created: 2025-12-22
% Last Modified: 2025-12-22
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
fs          = 5;                % 采样频率 (Hz)
t_total     = 2000;             % 总仿真时间长度 (秒)
num_samples = fs * t_total;     % 样本点数
time_axis   = (0:num_samples-1)' / fs;

% 输入信号参数（弱周期信号）
f0      = 0.01;               % 目标信号频率 (Hz)，满足 f0 << fs
A0      = 0.05;               % 信号幅值
clean_signal = A0 * sin(2*pi*f0*time_axis);

% 噪声强度 D 扫描范围
D_list   = 0.05:0.01:0.45;
num_D    = length(D_list);

% 每个 D 重复次数（用于统计平均）
n_repeat = 10;

% 预分配结果存储
snr_mean    = zeros(num_D, 1);
mpe_mean    = zeros(num_D, 1);
ami_mean    = zeros(num_D, 1);

%% 2. 定义双稳系统漂移函数 ==========================================
% drift_func = @(x) CBSR_Dynamics(x, 1, 1);
% drift_func = @(x) UBSR_Dynamics(x, 1, 1);
[a, b, k1, k2] = CalibrateHSUBSR(0.9, 3.6, 35);
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);

%% 3. 扫描噪声强度 D，统计 SNR 与 RSCM 指标 ==========================
fprintf('开始扫描噪声强度 D，共 %d 个取值，每个重复 %d 次...\n', num_D, n_repeat);
for idx_D = 1:num_D
    D = D_list(idx_D);
    fprintf('  -> processing D = %.4f ...\n', D);
    
    snr_rep     = zeros(n_repeat, 1);
    mpe_rep     = zeros(n_repeat, 1);
    ami_rep     = zeros(n_repeat, 1);
    
    parfor k = 1:n_repeat
        % ---- 生成噪声序列 (Gaussian white noise) ---------------------
        noise_seq = sqrt(2*D*fs) * randn(num_samples, 1);  % = sqrt(2D*fs)*randn
        % ---- 求解 SDE 获得系统输出 x(t) ------------------------------
        x = RK4Solver2(drift_func, clean_signal, noise_seq, fs);
        
        % ---- 去掉前 10% 瞬态，保留稳态部分 ---------------------------
        steady_start = round(0.1 * num_samples);
        x_stead = x(steady_start+1:end);
        
        % ---- 计算输出 SNR 指标 ---------------------
        snr_rep(k) = SNRo2(x_stead, fs, f0);
        
        % ---- 计算改进 AMI 指标（驻留时间自适应裁剪） -------
        [ami_rep(k), ~] = AMI2(x_stead, fs);

        % ---- 计算多尺度排列熵 MPE 指标 -----------------------------
        [scales, info] = SelectMpeScales(x_stead, fs, 3);
        [~, mpnorm, ~, ~] = MultiScalePermEn(x_stead, scales);
        mpe_rep(k) = std(mpnorm(~isnan(mpnorm)));  % 用 MPE 标准差作为指标    
    end
    
    % ---- 对重复试验求平均与标准差 -----------------------------------
    snr_mean(idx_D)     = mean(snr_rep);
    mpe_mean(idx_D)     = mean(mpe_rep);
    ami_mean(idx_D)     = mean(ami_rep);
end

%% 4. 计算峰值位置与相关系数 =========================================

% ---- 峰值噪声强度 D_peak --------------------------------------------
[~, idx_peak_snr]  = max(snr_mean);
[~, idx_peak_ami]  = max(ami_mean);
[~, idx_peak_mpe]  = max(1 - mpe_mean);

D_peak_snr      = D_list(idx_peak_snr);
D_peak_ami      = D_list(idx_peak_ami);
D_peak_mpe      = D_list(idx_peak_mpe);

delta_D_ami = abs(D_peak_ami - D_peak_snr);
delta_D_mpe = abs(D_peak_mpe - D_peak_snr);

% ---- 曲线相关系数 corr(SNR, AMI) -----------------------------------
R_ami = corrcoef(snr_mean, ami_mean);
rho_ami = R_ami(1, 2);

% ---- 曲线相关系数 corr(SNR, MPE) -----------------------------------
R_mpe = corrcoef(snr_mean, 1 - mpe_mean);
rho_mpe = R_mpe(1, 2);

fprintf('\n================ 指标对比结果 ================\n');
fprintf('SNR : 峰值 D_peak = %.4f\n', D_peak_snr);
fprintf('AMI : 峰值 D_peak = %.4f, |ΔD| = %.4f, corr(SNR, AMI) = %.3f\n', ...
    D_peak_ami, delta_D_ami, rho_ami);
fprintf('MPE : 峰值 D_peak = %.4f, |ΔD| = %.4f, corr(SNR, MPE) = %.3f\n', ...
    D_peak_mpe, delta_D_mpe, rho_mpe);
fprintf('===============================================\n');

%% 5. 绘制原始曲线对比 ================================================
SetThesisDefaultStyle();
fig = CreateThesisFigure();
subplot(2, 2, 1);
plot(D_list, smooth(snr_mean, 1), 'o-', 'DisplayName', 'SNR');
xlabel('D');
ylabel('SNR');

subplot(2, 2, 2);
plot(D_list, smooth(ami_mean, 1), 'd-.', 'DisplayName', 'AMI');
xlabel('D');
ylabel('AMI');

subplot(2, 2, 3);
plot(D_list, smooth(1 - mpe_mean, 1), 's--', 'DisplayName', 'MPE');
xlabel('D');
ylabel('1-MPE');

subplot(2, 2, 4);
plot(D_list, smooth((1 - mpe_mean).*ami_mean, 1), 'd-.', 'DisplayName', 'RSCM'); 
xlabel('D');
ylabel('RSCM');

results.snr = snr_mean;
results.ami = ami_mean;
results.mpe = mpe_mean;
results.D_list = D_list;