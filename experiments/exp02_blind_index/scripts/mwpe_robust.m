% =========================================================================
% Description: 验证自适应尺度选择策略超参数对MWPE指标的影响
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

mpe_list2 = zeros(n_D, 1); % 存储每次MC的MPE结果
mpe_list4 = zeros(n_D, 1);
mpe_list6 = zeros(n_D, 1);

snr_list = zeros(n_D, 1);

h = 1 /fs;
%% 2. 循环计算
m = 3;  % MPE嵌入维度
parfor i = 1 : n_D
    D = D_list(i);
    mpe_rep2 = zeros(num_mc, 1);
    mpe_rep4 = zeros(num_mc, 1);
    mpe_rep6 = zeros(num_mc, 1);
    snr_rep = zeros(num_mc, 1);
    
    for mc_iter = 1:num_mc
        % 1. 生成含噪信号
        noise = sqrt(2*D*fs) * randn(N, 1);
        input = clean_sig + noise;
        
        % 2. RK4 求解 SR 输出
        x = RK4Solver2(drift_func, clean_sig, noise, fs);
        
        % 3. 预处理 (去瞬态 + 归一化)
        x_steady = x(round(0.1*N)+1:end);
        
        out2 = AdaptiveScalesACF(x_steady, fs, 'm', m, 'C', 50);
        [~, mpe2, ~, ~]  = MultiScalePermEn(x_steady, out2.S, 'm', m);
        H_min2 = mean(mpe2);
        
        % out4 = AdaptiveScalesACF(x_steady, fs, 'm', m, 'C', 4);
        % [~, mpe4, ~, ~]  = MultiScalePermEn(x_steady, out4.S, 'm', m);
        % H_min4 = mean(mpe4);
        
        % out6 = AdaptiveScalesACF(x_steady, fs, 'm', m, 'C', 6);
        % [~, mpe6, ~, ~]  = MultiScalePermEn(x_steady, out6.S, 'm', m);
        % H_min6 = mean(mpe6);
        
        snr_rep(mc_iter) = SNRo2(x_steady, fs, f0);
        
        % 累积求均值/方差
        mpe_rep2(mc_iter) = H_min2;
        % mpe_rep4(mc_iter) = H_min4;
        % mpe_rep6(mc_iter) = H_min6;
    end
    
    mpe_list2(i) = mean(mpe_rep2);
    % mpe_list4(i) = mean(mpe_rep4);
    % mpe_list6(i) = mean(mpe_rep6);
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
plot(D_list, smooth(mpe_list2, 1), '-^', 'LineWidth', 2); hold on;
% plot(D_list, smooth(mpe_list4, 1), '-s', 'LineWidth', 2);
% plot(D_list, smooth(mpe_list6, 1), '-d', 'LineWidth', 2);

results.input.fs = fs;
results.input.T = T;
results.input.A0 = A0;
results.input.f0 = f0;
results.input.clean_sig = clean_sig;
results.potential.xm = 1;
results.potential.dU = 0.25;
results.potential.shape = 1.01;
results.output.snr_list = snr_list;
results.output.mpe_list2 = mpe_list2;
results.output.mpe_list4 = mpe_list4;
results.output.mpe_list6 = mpe_list6;