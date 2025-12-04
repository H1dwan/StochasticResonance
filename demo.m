% =========================================================================
% Description: 代码库使用示例脚本，验证 CBSR 的输出信噪比曲线随噪声强度间的关系
%
% Author: LiuShuang
% Created: 2025-12-04
% Last Modified: 2025-12-04
%
% Usage: Just run this script
%
% Input:
%   None
% Output:
%   SNR vs D 曲线
% =========================================================================

clc; clear; close all;

%% 1. 基本仿真参数设置 ========================================

fs  = 20;          % 采样频率 Hz
T   = 1000;        % 信号时长 s
N   = fs * T;      % 样本点数
t   = (0:N-1)'/fs; % 时间向量

n_repeat = 10;     % 重复仿真次数（多次平均减小随机性）

A0 = 0.05;         % 弱信号幅值
f0 = 0.01;         % 弱信号频率 Hz
s = A0 * sin(2*pi*f0*t);

a = 1;            % CBSR 参数 a
b = 1;             % CBSR 参数 b

xm      = sqrt(a / b);      % 势阱位置 ±xm
DeltaU  = a^2 / (4 * b);    % 势垒高度 ΔU_c = a^2/(4b)
drift   = @(x) CBSR_Dynamics(x, a, b);

fprintf('CBSR: xm = %.4f, DeltaU = %.4f\n', xm_c, DeltaU_c);

%% 2. 不同噪声强度 D 下的 SNR 曲线（A 固定） ===================
D_list = 0.01:0.01:1;
snr_list = zeros(length(D_list), 1);

for iD = 1:length(D_list)
    D = D_list(iD);
    snr_rep = zeros(nReal,1);
    fprintf("processing D = %.3f ...\n", D);
    
    for k = 1:n_repeat
        % ---- 生成噪声序列 noise_seq ----------------
        h = 1/fs;
        noise_seq = sqrt(2*D/h) * randn(N,1);   % = sqrt(2D*fs)*randn
        
        % ---- 求解 HESUBSR 系统输出 -----------------
        x = RK4Solver2(drift_CBSR, s, noise_seq, fs);
        
        % ---- 去掉前一段瞬态，避免初始条件影响 ----
        steady_start = round(0.2*N);
        x_stead = x(steady_start:end);
        
        % ---- 计算输出 SNR -------------------------
        snr_rep(k) = SNRo(x_stead, fs, f0);
    end
    
    % 多次仿真求平均 SNR
    snr_list(iD) = mean(snr_rep);
end

%% 3. 绘制 SNR - D 曲线 ========================================
figure('Position', [100, 100, 1200, 900], 'Color', 'white');
plot(D_list, snr_list, 'bo-', 'LineWidth', 1.5);
xlabel('噪声强度 D');
ylabel('信噪比 SNR (dB)');
title(sprintf('SNR-D 曲线 (A = %.3f)', A0));