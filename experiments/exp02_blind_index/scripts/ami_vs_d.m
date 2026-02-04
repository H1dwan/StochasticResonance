% =========================================================================
% Description: 探究 AMI 指标在不同噪声强度下的演化过程
%
% Author: LiuShuang
% Created: 2026-01-24
% Last Modified: 2026-01-24
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 仿真参数设置
fs  = 5;              % 采样率
T   = 2000;           % 信号时长
N   = fs * T;
t   = (0:N-1)'/fs;

% 弱信号参数
A0  = 0.1;
f0  = 0.01;
clean_sig = A0 * sin(2*pi*f0*t);

% SR 系统参数 (a=1, b=1)
[a, b, k1, k2] = CalibrateHSUBSR(1, 0.25, 1.01);
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
% drift_func = @(x) CBSR_Dynamics(x, 1, 1);

% 扫描范围
D_list  = 0.05:0.01:0.45;   % 噪声强度范围
% D_list = 0.12;
n_D     = length(D_list);
num_mc  = 100;              % 蒙特卡洛次数，可按需调整

ami_list = zeros(n_D, 1);   % 存储每次MC的AMI结果

snr_list = zeros(n_D, 1);

h = 1 /fs;

%% 2. 循环计算
parfor i = 1 : n_D
    D = D_list(i);
    ami_rep = zeros(num_mc, 1);
    snr_rep = zeros(num_mc, 1);
    
    for mc_iter = 1:num_mc
        % 1. 生成含噪信号
        noise = sqrt(2*D*fs) * randn(N, 1);
        input = clean_sig + noise;
        
        % 2. RK4 求解 SR 输出
        x = RK4Solver2(drift_func, clean_sig, noise, fs);
        
        % 3. 预处理 (去瞬态 + 归一化)
        x_steady = x(round(0.1*N)+1:end);
        
        [ami_rep(mc_iter), ~] = AMI2(x_steady, fs);
        
        snr_rep(mc_iter) = SNRo2(x_steady, fs, f0);
        
    end
    
    ami_list(i) = mean(ami_rep);

    snr_list(i) = mean(snr_rep);
    
    if mod(i, 5) == 0
        fprintf('  D index %d/%d 完成\n', i, n_D);
    end
end

%% 绘制
figure;
yyaxis left;
plot(D_list, smooth(snr_list, 1), '-o', 'LineWidth', 2);
yyaxis right;
plot(D_list, smooth(ami_list, 1), '-^', 'LineWidth', 2); hold on;


results.input.fs = fs;
results.input.T = T;
results.input.A0 = A0;
results.input.f0 = f0;
results.input.clean_sig = clean_sig;
results.potential.xm = 1;
results.potential.dU = 0.25;
results.potential.shape = 1.01;
results.output.snr_list = snr_list;
results.output.ami_list = ami_list;
