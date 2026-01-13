% =========================================================================
% Description: 公平对比实验 (Fair Comparison)
%              统一采用 (xm, dU) 作为几何特征参数进行 PSO 优化
%
% Models:
%   1. UBSR   (xm, dU) -> (a, b)
%   2. PLBSR  (xm, dU) -> (L0, U0)
%   3. HSUBSR (xm, dU, shape) -> (a, b, k1, k2)
%
% Author: Stochastic Resonance Research Assistant
% =========================================================================

clc; clear; close all;

%% 1. 信号与噪声环境设置
% -------------------------------------------------------------------------
rng(42);                % 固定随机种子
fs  = 5;                % 采样频率 Hz
T   = 1000;             % 信号时长 s
N   = fs * T;           % 采样点数
t   = (0:N-1)'/fs;

A0  = 0.1;              % 弱信号幅值
f0  = 0.01;             % 信号频率 Hz
D   = 0.6;              % 噪声强度

% 构造输入信号
pure_sig  = A0 * sin(2*pi*f0*t);
noise_seq = sqrt(2*D*fs) * randn(N, 1);
input_sig = pure_sig + noise_seq;
snr_in = SNRo2(input_sig, fs, f0);

fprintf('=== 公平对比仿真开始 (基于几何参数 xm, dU) ===\n');
fprintf('信号: A=%.2f, D=%.2f, SNR=%.2f dB\n', A0, D, snr_in);

%% 2. 统一优化参数设置
% -------------------------------------------------------------------------
pop_size = 10;
max_iter = 20;

% --- 几何参数的搜索边界 (对所有模型统一) ---
% xm (势阱位置): [0.1, 3.0]
% dU (势垒高度): [0.1, 3.0]
xm_bound = [0.1, 3.0];
dU_bound = [0.1, 3.0];

% 1. UBSR 参数: [xm, dU]
lb_ubsr = [xm_bound(1), dU_bound(1)];
ub_ubsr = [xm_bound(2), dU_bound(2)];
dim_ubsr = 2;

% 2. PLBSR 参数: [xm, dU]
lb_plbsr = [xm_bound(1), dU_bound(1)];
ub_plbsr = [xm_bound(2), dU_bound(2)];
dim_plbsr = 2;

% 3. HSUBSR 参数: [xm, dU, shape]
% shape 因子额外设定范围
lb_hs = [xm_bound(1), dU_bound(1), 1.1]; 
ub_hs = [xm_bound(2), dU_bound(2), 50];
dim_hs = 3;

%% 3. 执行 PSO 优化
% -------------------------------------------------------------------------
% 目标函数句柄封装
fobj_ubsr  = @(p) Fitness_UBSR_Fair(p, pure_sig, noise_seq, fs, f0);
fobj_plbsr = @(p) Fitness_PLBSR_Fair(p, pure_sig, noise_seq, fs, f0);
fobj_hs    = @(p) Fitness_HSUBSR_Fair(p, pure_sig, noise_seq, fs, f0);

fprintf('\n(1/3) 正在优化 UBSR ...\n');
[best_score_u, best_pos_u, curve_u] = PSO(pop_size, max_iter, lb_ubsr, ub_ubsr, dim_ubsr, fobj_ubsr, false);
best_snr_u = -best_score_u;
[a_opt, b_opt] = MapParams_UBSR(best_pos_u(1), best_pos_u(2));
fprintf('   -> 最佳 SNR: %.4f dB | 几何参数: xm=%.4f, dU=%.4f | 物理参数: a=%.4f, b=%.4f\n', ...
    best_snr_u, best_pos_u(1), best_pos_u(2), a_opt, b_opt);

fprintf('\n(2/3) 正在优化 PLBSR ...\n');
[best_score_p, best_pos_p, curve_p] = PSO(pop_size, max_iter, lb_plbsr, ub_plbsr, dim_plbsr, fobj_plbsr, false);
best_snr_p = -best_score_p;
fprintf('   -> 最佳 SNR: %.4f dB | 几何参数: xm=%.4f, dU=%.4f\n', ...
    best_snr_p, best_pos_p(1), best_pos_p(2));

fprintf('\n(3/3) 正在优化 HSUBSR ...\n');
[best_score_h, best_pos_h, curve_h] = PSO(pop_size, max_iter, lb_hs, ub_hs, dim_hs, fobj_hs, false);
best_snr_h = -best_score_h;
fprintf('   -> 最佳 SNR: %.4f dB | 几何参数: xm=%.4f, dU=%.4f, Shape=%.1f\n', ...
    best_snr_h, best_pos_h(1), best_pos_h(2), best_pos_h(3));

%% 4. 结果可视化
% -------------------------------------------------------------------------
% 计算最佳波形
x_ubsr  = GetOutput_UBSR_Fair(best_pos_u, pure_sig, noise_seq, fs);
x_plbsr = GetOutput_PLBSR_Fair(best_pos_p, pure_sig, noise_seq, fs);
x_hs    = GetOutput_HSUBSR_Fair(best_pos_h, pure_sig, noise_seq, fs);

figure('Position', [100, 100, 1000, 500], 'Color', 'white');

% 子图1：收敛曲线
subplot(1, 2, 1);
plot(1:max_iter, -curve_u, 'g--', 'LineWidth', 1.5); hold on;
plot(1:max_iter, -curve_p, 'b-.', 'LineWidth', 1.5);
plot(1:max_iter, -curve_h, 'r-', 'LineWidth', 2.0);
xlabel('Iter'); ylabel('Max SNR (dB)');
title('PSO 优化收敛曲线对比');
legend('UBSR', 'PLBSR', 'HSUBSR', 'Location', 'SouthEast');
grid on;

% 子图2：时域波形片段
subplot(1, 2, 2);
idx = 2000:3000; % 选取中间一段
plot(t(idx), x_ubsr(idx), 'g--', 'LineWidth', 1); hold on;
plot(t(idx), x_plbsr(idx), 'b-.', 'LineWidth', 1);
plot(t(idx), x_hs(idx), 'r-', 'LineWidth', 1.5);
title('最优输出波形对比');
xlabel('Time (s)');
legend('UBSR', 'PLBSR', 'HSUBSR');
grid on;

% 绘制输出时频图
Plot_Time_Frequency(x_ubsr, fs, length(x_ubsr));
Plot_Time_Frequency(x_plbsr, fs, length(x_plbsr));
Plot_Time_Frequency(x_hs, fs, length(x_hs));

%% 6. 保存结果
results.input.A0 = A0;
results.input.D  = D;
results.input.T = T;
results.input.fs = fs;
results.input.f0 = f0;
results.input.pure_sig = pure_sig;
results.input.noise_seq = noise_seq;
results.best_params.UBSR = best_pos_u;
results.best_params.PLBSR = best_pos_p;
results.best_params.HSUBSR = best_pos_h;
results.outputs.UBSR = x_ubsr;
results.outputs.PLBSR = x_plbsr;
results.outputs.HSUBSR = x_hs;

%% ---------------- 核心辅助函数 (参数映射) ---------------- %%

% --- 1. UBSR 映射与适应度 ---
function [a, b] = MapParams_UBSR(xm, dU)
    % 物理推导: U(x) = -a/2 x^2 + b/4 x^4
    % xm^2 = a/b; dU = a^2/(4b)
    % 解得: a = 4*dU / xm^2; b = 4*dU / xm^4
    a = 4 * dU / (xm^2);
    b = 4 * dU / (xm^4);
end

function val = Fitness_UBSR_Fair(p, s, n, fs, f0)
    [a, b] = MapParams_UBSR(p(1), p(2)); % p=[xm, dU]
    drift = @(x) UBSR_Dynamics(x, a, b);
    x = RK4Solver2(drift, s, n, fs);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, fs, f0);
end

function x = GetOutput_UBSR_Fair(p, s, n, fs)
    [a, b] = MapParams_UBSR(p(1), p(2));
    drift = @(x) UBSR_Dynamics(x, a, b);
    x = RK4Solver2(drift, s, n, fs);
end

% --- 2. PLBSR 映射与适应度 ---
function [U0, L0] = MapParams_PLBSR(xm, dU)
    % PLBSR 参数直接对应: L0是阱宽(xm), U0是垒高(dU)
    L0 = xm;
    U0 = dU;
end

function val = Fitness_PLBSR_Fair(p, s, n, fs, f0)
    [U0, L0] = MapParams_PLBSR(p(1), p(2)); % p=[xm, dU]
    drift = @(x) PLBSR_Dynamics(x, U0, L0);
    x = RK4Solver2(drift, s, n, fs);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, fs, f0);
end

function x = GetOutput_PLBSR_Fair(p, s, n, fs)
    [U0, L0] = MapParams_PLBSR(p(1), p(2));
    drift = @(x) PLBSR_Dynamics(x, U0, L0);
    x = RK4Solver2(drift, s, n, fs);
end

% --- 3. HSUBSR 映射与适应度 ---
function val = Fitness_HSUBSR_Fair(p, s, n, fs, f0)
    % p = [xm, dU, shape]
    % 调用已有的校准函数
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    x = RK4Solver2(drift, s, n, fs);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, fs, f0);
end

function x = GetOutput_HSUBSR_Fair(p, s, n, fs)
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    x = RK4Solver2(drift, s, n, fs);
end