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
% rng(1);            % 固定随机种子，保证可重复
fs  = 5;           % 采样频率 Hz
T   = 2000;        % 信号时长 s
N   = fs * T;      % 样本点数
t   = (0:N-1)'/fs; % 时间向量

n_repeat = 10;     % 重复仿真次数（多次平均减小随机性）

A0 = 0.05;         % 弱信号幅值
f0 = 0.01;         % 弱信号频率 Hz
s = A0 * sin(2*pi*f0*t);

a = 1;            % CBSR 参数 a
b = 1;             % CBSR 参数 b

xm      = sqrt(a / b);      % 势阱位置 ±xm
dU  = a^2 / (4 * b);    % 势垒高度 ΔU_c = a^2/(4b)

drift   = @(x) CBSR_Dynamics(x, a, b);
% drift = @(x) HESUBSR_Dynamics(x, 1.1136, 1, 1, 1.6703);
% shape = 20;
% [dslc.a, dslc.b, dslc.k1, dslc.k2] = CalibrateHSUBSR(xm, dU, shape);
% drift = @(x) HSUBSR_Dynamics(x, dslc.a, dslc.b, dslc.k1, dslc.k2);

fprintf('CBSR: xm = %.4f, DeltaU = %.4f\n', xm, dU);

%% 2. 不同噪声强度 D 下的 SNR 曲线（A 固定） ===================
D_list = 0.05:0.01:0.45;
snr_list = zeros(length(D_list), 1);
metric1_curve = zeros(length(D_list), 1);
metric2_curve = zeros(length(D_list), 1);

for iD = 1:length(D_list)
    D = D_list(iD);
    snr_rep = zeros(n_repeat,1);
    metric1_rep = zeros(n_repeat,1);
    metric2_rep = zeros(n_repeat,1);
    fprintf("processing D = %.3f ...\n", D);
    
    parfor k = 1:n_repeat
        % ---- 生成噪声序列 noise_seq ----------------
        noise_seq = sqrt(2*D*fs) * randn(N,1);   % = sqrt(2D*fs)*randn
        
        % ---- 求解 HESUBSR 系统输出 -----------------
        x = RK4Solver2(drift, s, noise_seq, fs);
        
        % ---- 去掉前一段瞬态，避免初始条件影响 ----
        steady_start = round(0.1*N);
        x_stead = x(steady_start+1:end);
        
        % ---- 计算输出 SNR -------------------------
        snr_rep(k) = SNRo2(x_stead, fs, f0);
        
        % ---- 计算新的指标（PermEn 可能返回向量，这里取第一个量化值） ----
        % [~, metric1_tmp, info1] = PermEn(x_stead, 'm', 4, 'tau', 75);
        % [~, metric2_tmp, info2] = PermEn(x_stead, 'm', 4, 'tau', 100);
        [metric1_rep(k), ~] = AMI2(x_stead, fs);
        % metric1_rep(k) = metric1_tmp(4);
        % metric2_rep(k) = metric2_tmp(4);
        % if isstruct(info1) && isfield(info1, 'f_star')
        %     fprintf("f star: %.4f Hz\n", info1.f_star);
        % end
    end
    
    % 多次仿真求平均 SNR
    snr_list(iD) = mean(snr_rep);
    metric1_curve(iD) = mean(metric1_rep);
    % metric2_curve(iD) = mean(metric2_rep);
end

%% 3. 绘制 SNR - D 曲线 ========================================
figure('Position', [100, 100, 1200, 900], 'Color', 'white');
plot(D_list, snr_list, 'bo-', 'LineWidth', 1.5); hold on;
xlabel('噪声强度 D');
ylabel('信噪比 SNR (dB)');
title(sprintf('SNR-D 曲线 (A = %.3f)', A0));

figure('Position', [100, 100, 1200, 900], 'Color', 'white');
plot(D_list, metric1_curve, 'ro-', 'LineWidth', 1.5, 'DisplayName', '75'); hold on;
% plot(D_list, metric2_curve, 'go-', 'LineWidth', 1.5, 'DisplayName', '100');
xlabel('噪声强度 D');
ylabel('归一化谱熵 SpecEntropy');
title(sprintf('谱熵-D 曲线 (A = %.3f)', A0));
legend('show');

SetThesisDefaultStyle();
CreateThesisFigure();
tiledlayout(1 , 1 ,'Padding','compact','TileSpacing','compact');
yyaxis left;
plot(D_list, snr_list, 'o-', 'LineWidth', 2);
ylabel('SNR')
xlabel('$D$')
yyaxis right;
plot(D_list, metric1_curve, 's--', 'LineWidth', 2);
ylabel('AMI')
legend('SNR', 'AMI');