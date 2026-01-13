% =========================================================================
% Description: Alpha稳定噪声下三种抗饱和模型(UBSR/PLBSR/HSUBSR)性能对比
%
% Author: LiuShuang
% Created: 2026-01-13
% Last Modified: 2026-01-13
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 仿真环境设置 =========================================================
fprintf('=== 仿真开始: Alpha 稳定噪声下的模型性能对比 ===\n');

% 1.1 基础参数
fs         = 5;               % 采样频率 Hz (绝热近似)
T          = 1000;            % 信号时长 s (足够长以保证统计特性)
n_samples  = fs * T;
t_vec      = (0:n_samples-1)' / fs;
h          = 1/fs;

% 1.2 弱信号参数
A0         = 0.1;             % 弱信号幅值
f0         = 0.01;            % 信号频率
clean_sig  = A0 * sin(2*pi*f0*t_vec);

% 1.3 Alpha 稳定噪声参数
alpha_noise = 1.2;            % 特征指数 < 2 (表现出显著脉冲性/厚尾)
beta_noise  = 0;              % 对称分布
sigma_noise = 2;            % 噪声强度 scale parameter
mu_noise    = 0;

fprintf('信号参数: A=%.2f, f=%.2f Hz\n', A0, f0);
fprintf('噪声参数: alpha=%.1f (脉冲噪声), sigma=%.1f\n', alpha_noise, sigma_noise);

% 1.4 生成噪声
% 注意: sigma_discrete 是离散增量的尺度参数
sigma_discrete = sigma_noise * h^(1/alpha_noise);
levy_increments = GenerateAlphaStableNoise(alpha_noise, beta_noise, sigma_discrete, mu_noise, n_samples, 1);

% 关键修正: RK4Solver2 会执行 dx += noise * h
% 因此传入的 noise_seq 应该是 "Force" = Increment / h
noise_force_for_solver = levy_increments / 1;

% 混合信号 (用于计算输入信噪比和绘图)
s_noisy = clean_sig + levy_increments;
snr_in = SNRo2(s_noisy, fs, f0);
fprintf('输入混合信号 SNR = %.4f dB\n', snr_in);

Plot_Time_Frequency(s_noisy, fs, length(s_noisy));

%% 2. 优化参数配置 (统一几何参数) ===========================================
% 统一搜索空间: 势阱位置 xm, 势垒高度 dU
% HSUBSR 额外搜索形状因子 S
pop_size = 30;
max_iter = 10;

% 搜索边界 [xm_min, xm_max]
bound_xm = [0.3, 5.0];
% 搜索边界 [dU_min, dU_max]
bound_dU = [0.1, 8.0];

% --- 模型 1: UBSR ---
lb_u = [bound_xm(1), bound_dU(1)];
ub_u = [bound_xm(2), bound_dU(2)];
dim_u = 2;
fobj_u = @(p) Fitness_UBSR_Fair(p, clean_sig, noise_force_for_solver, fs, f0);

% --- 模型 2: PLBSR ---
lb_p = [bound_xm(1), bound_dU(1)];
ub_p = [bound_xm(2), bound_dU(2)];
dim_p = 2;
fobj_p = @(p) Fitness_PLBSR_Fair(p, clean_sig, noise_force_for_solver, fs, f0);

% --- 模型 3: HSUBSR ---
% 第三个参数为形状因子 S (1.1~50)
lb_h = [bound_xm(1), bound_dU(1), 1.1];
ub_h = [bound_xm(2), bound_dU(2), 100];
dim_h = 3;
fobj_h = @(p) Fitness_HSUBSR_Fair(p, clean_sig, noise_force_for_solver, fs, f0);

%% 3. 执行 PSO 寻优 ========================================================
fprintf('\n正在优化 UBSR 模型...\n');
[best_score_u, best_pos_u, curve_u] = PSO(pop_size, max_iter, lb_u, ub_u, dim_u, fobj_u, false);
best_snr_u = -best_score_u;

fprintf('正在优化 PLBSR 模型...\n');
[best_score_p, best_pos_p, curve_p] = PSO(pop_size, max_iter, lb_p, ub_p, dim_p, fobj_p, false);
best_snr_p = -best_score_p;

fprintf('正在优化 HSUBSR 模型...\n');
[best_score_h, best_pos_h, curve_h] = PSO(pop_size, max_iter, lb_h, ub_h, dim_h, fobj_h, false);
best_snr_h = -best_score_h;

%% 4. 结果计算与展示 =======================================================
% 重新计算最优输出
x_ubsr  = GetOutput_UBSR(best_pos_u, clean_sig, noise_force_for_solver, fs);
x_plbsr = GetOutput_PLBSR(best_pos_p, clean_sig, noise_force_for_solver, fs);
x_hsubsr= GetOutput_HSUBSR(best_pos_h, clean_sig, noise_force_for_solver, fs);

% 打印对比表
fprintf('\n======================================================\n');
fprintf('Alpha 稳定噪声 (alpha=%.1f) 下的性能对比\n', alpha_noise);
fprintf('======================================================\n');
fprintf('| 模型    | 最佳 SNR (dB) | 提升 (vs UBSR) | 最优参数 (xm, dU)\n');
fprintf('|---------|---------------|----------------|-------------------\n');
fprintf('| UBSR    | %10.4f    |       -        | (%.2f, %.2f)\n', best_snr_u, best_pos_u(1), best_pos_u(2));
fprintf('| PLBSR   | %10.4f    | %+10.4f dB   | (%.2f, %.2f)\n', best_snr_p, best_snr_p-best_snr_u, best_pos_p(1), best_pos_p(2));
fprintf('| HSUBSR  | %10.4f    | %+10.4f dB   | (%.2f, %.2f, S=%.1f)\n', best_snr_h, best_snr_h-best_snr_u, best_pos_h(1), best_pos_h(2), best_pos_h(3));
fprintf('======================================================\n');

% --- 绘图 1: 收敛曲线 ---
figure('Position', [100, 100, 1000, 400], 'Color', 'white');
subplot(1, 2, 1);
plot(1:max_iter, -curve_u, 'g--', 'LineWidth', 1.5); hold on;
plot(1:max_iter, -curve_p, 'b-.', 'LineWidth', 1.5);
plot(1:max_iter, -curve_h, 'r-', 'LineWidth', 2.0);
xlabel('迭代次数'); ylabel('最佳 SNR (dB)');
title(['PSO 收敛曲线 (\alpha = ' num2str(alpha_noise) ')']);
legend('UBSR', 'PLBSR', 'HSUBSR');
grid on;

% --- 绘图 2: 时域波形 (局部放大) ---
subplot(1, 2, 2);
range_idx = 2000:3000; % 取中间一段
t_sub = t_vec(range_idx);
plot(t_sub, x_ubsr(range_idx), 'g--', 'LineWidth', 1); hold on;
plot(t_sub, x_plbsr(range_idx), 'b-.', 'LineWidth', 1);
plot(t_sub, x_hsubsr(range_idx), 'r-', 'LineWidth', 1.5);
% 绘制等效的脉冲噪声位置 (示意)
large_noise_idx = find(abs(levy_increments(range_idx)) > 3*std(levy_increments));
if ~isempty(large_noise_idx)
    scatter(t_sub(large_noise_idx), zeros(size(large_noise_idx)), 'kx', 'DisplayName', 'Large Noise Impulses');
end
xlabel('时间 (s)'); ylabel('输出幅值');
title('最优输出波形片段');
legend('UBSR', 'PLBSR', 'HSUBSR');
grid on;

Plot_Time_Frequency(x_ubsr, fs, length(x_ubsr));
Plot_Time_Frequency(x_plbsr, fs, length(x_plbsr));
Plot_Time_Frequency(x_hsubsr, fs, length(x_hsubsr));

%% 5. 保存结果 ============================================================
results.inputs.fs = fs;
results.inputs.T = T;
results.inputs.A0 = A0;
results.inputs.f0 = f0;
results.inputs.clean_sig = clean_sig;
results.inputs.noise_force = noise_force_for_solver;
results.outputs.UBSR = x_ubsr;
results.outputs.PLBSR = x_plbsr;
results.outputs.HSUBSR = x_hsubsr;
results.best_params.UBSR = best_pos_u;
results.best_params.PLBSR = best_pos_p;
results.best_params.HSUBSR = best_pos_h;

%% 辅助函数定义 ============================================================

% --- 1. 参数映射与适应度函数 (UBSR) ---
function [a, b] = Map_UBSR(xm, dU)
a = 4 * dU / (xm^2);
b = 4 * dU / (xm^4);
end
function val = Fitness_UBSR_Fair(p, s, n, fs, f0)
[a, b] = Map_UBSR(p(1), p(2));
drift = @(x) UBSR_Dynamics(x, a, b);
x = RK4Solver2(drift, s, n, fs);
% 注意: Alpha噪声可能导致发散，加一个保护
if any(isnan(x)) || any(isinf(x)), val = 1e5; return; end
val = -SNRo2(x(round(0.1*end):end), fs, f0);
end
function x = GetOutput_UBSR(p, s, n, fs)
[a, b] = Map_UBSR(p(1), p(2));
drift = @(x) UBSR_Dynamics(x, a, b);
x = RK4Solver2(drift, s, n, fs);
end

% --- 2. 参数映射与适应度函数 (PLBSR) ---
function [U0, L0] = Map_PLBSR(xm, dU)
L0 = xm; U0 = dU;
end
function val = Fitness_PLBSR_Fair(p, s, n, fs, f0)
[U0, L0] = Map_PLBSR(p(1), p(2));
drift = @(x) PLBSR_Dynamics(x, U0, L0);
x = RK4Solver2(drift, s, n, fs);
if any(isnan(x)) || any(isinf(x)), val = 1e5; return; end
val = -SNRo2(x(round(0.1*end):end), fs, f0);
end
function x = GetOutput_PLBSR(p, s, n, fs)
[U0, L0] = Map_PLBSR(p(1), p(2));
drift = @(x) PLBSR_Dynamics(x, U0, L0);
x = RK4Solver2(drift, s, n, fs);
end

% --- 3. 参数映射与适应度函数 (HSUBSR) ---
function val = Fitness_HSUBSR_Fair(p, s, n, fs, f0)
[a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x = RK4Solver2(drift, s, n, fs);
if any(isnan(x)) || any(isinf(x)), val = 1e5; return; end
val = -SNRo2(x(round(0.1*end):end), fs, f0);
end
function x = GetOutput_HSUBSR(p, s, n, fs)
[a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x = RK4Solver2(drift, s, n, fs);
end

% --- 4. Alpha 稳定噪声生成器 (CMS Algorithm) ---
function x = GenerateAlphaStableNoise(alpha, beta, sigma, mu, m, n)
if nargin < 6, n = 1; end
if nargin < 5, m = 1; end
if nargin < 4, mu = 0; end
if nargin < 3, sigma = 1; end

% CMS 算法实现
V = (rand(m, n) - 0.5) * pi;
W = exprnd(1, m, n);

if alpha == 2
    x = sigma * sqrt(2) * randn(m, n) + mu;
    return;
end
if alpha == 1
    term1 = (pi/2 + beta * V) .* tan(V);
    arg_log = (pi/2 * W .* cos(V)) ./ (pi/2 + beta * V);
    arg_log(arg_log <= 0) = eps;
    term2 = beta * log(arg_log);
    X = 2/pi * (term1 - term2);
    x = sigma * X + (2/pi * beta * sigma * log(sigma)) + mu;
    return;
end

B = atan(beta * tan(pi * alpha / 2)) / alpha;
S = (1 + (beta * tan(pi * alpha / 2))^2)^(1 / (2 * alpha));
num = sin(alpha * (V + B));
den = cos(V).^(1 / alpha);
term_inner = cos(V - alpha * (V + B)) ./ W;
term_inner(term_inner <= 0) = eps;
term2 = term_inner.^((1 - alpha) / alpha);
Z = S * (num ./ den) .* term2;
x = sigma * Z + mu;
end