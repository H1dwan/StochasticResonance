% =========================================================================
% Description: 三层级联随机共振系统在 Alpha 稳定噪声下的检测能力验证
%
% Author: LiuShuang
% Created: 2026-03-03
% Last Modified: 2026-03-03
% =========================================================================

clc; clear; close all;
% rng(0);  % 固定随机种子，便于复现

%% 1. 仿真参数设置 ========================================================
fs = 5;                  % 采样频率 Hz
T  = 1000;               % 信号时长 s
N  = fs * T;
t  = (0:N-1)' / fs;
h  = 1/fs;

% 弱周期信号
f0 = 0.01;               % 目标频率 Hz
A0 = 0.10;               % 信号幅值
clean_sig = A0 * sin(2*pi*f0*t);

% Alpha 稳定噪声参数
alpha_noise = 1.2;       % 特征指数 (<2 时为脉冲厚尾噪声)
beta_noise  = 0;         % 对称分布
sigma_noise = 2.0;       % 尺度参数
mu_noise    = 0;

% 级联配置
num_layers = 3;          % 三层级联

% Monte Carlo 配置
num_trials = 1;         % 重复试验次数

% PSO 配置
search_agents = 20;
max_iter = 50;
lb = [0.1, 0.05, 1.01];   % [xm, dU, shape]
ub = [2.0, 2.00, 1.99];

fprintf('=== 三层级联 SR 在 Alpha 噪声下检测能力验证 ===\n');
fprintf('fs=%d Hz, T=%d s, f0=%.3f Hz, A0=%.3f\n', fs, T, f0, A0);
fprintf('alpha=%.2f, sigma=%.2f, MonteCarlo=%d\n\n', alpha_noise, sigma_noise, num_trials);

%% 2. Monte Carlo 主循环 ==================================================
input_snr_list = zeros(num_trials, 1);
snr_1layer_list = zeros(num_trials, 1);
snr_3layer_list = zeros(num_trials, 1);

params_1layer = zeros(num_trials, 3);
params_3layer = zeros(num_trials, 3*num_layers);
input_noisy_all = zeros(N, num_trials);
layer_output_all = zeros(N, num_layers, num_trials);
layer_snr_all = zeros(num_trials, num_layers);

for trial = 1:num_trials
    % 2.1 生成 Alpha 稳定噪声增量
    sigma_discrete = sigma_noise * h^(1/alpha_noise);
    levy_inc = GenerateAlphaStableNoise(alpha_noise, beta_noise, sigma_discrete, mu_noise, N, 1);
    
    % RK4Solver2 中噪声项为: noise_seq(i) * h
    % 因此传入“随机力”应为增量/步长
    noise_force = levy_inc / 1;
    
    % 输入观测信号（用于计算输入检测指标）
    noisy_input = clean_sig + levy_inc;
    input_noisy_all(:, trial) = noisy_input;
    input_snr_list(trial) = SNRo2(noisy_input, fs, f0);
    
    % 2.2 单层 HSUBSR 基线（同噪声条件）
    fobj_1 = @(theta) CostFunc_HSUBSR_Layer(theta, clean_sig, noise_force, fs, f0);
    [best_val_1, best_pos_1, ~] = PSO(search_agents, max_iter, lb, ub, 3, fobj_1, false);
    params_1layer(trial, :) = best_pos_1;
    snr_1layer_list(trial) = -best_val_1;
    
    % 2.3 三层级联 HSUBSR（全局参数寻优）
    dim_global = 3 * num_layers;
    lb_global = repmat(lb, 1, num_layers);
    ub_global = repmat(ub, 1, num_layers);
    
    fobj_global = @(theta_vec) CostFunc_HSUBSR_Global(theta_vec, clean_sig, noise_force, fs, f0, num_layers);
    [best_val_global, best_pos_global, ~] = PSO(search_agents, max_iter, lb_global, ub_global, dim_global, fobj_global, false);
    
    [layer_outputs, layer_snrs] = SimulateHSUBSRCascade(best_pos_global, clean_sig, noise_force, fs, f0, num_layers);
    
    params_3layer(trial, :) = best_pos_global;
    layer_output_all(:, :, trial) = layer_outputs;
    layer_snr_all(trial, :) = layer_snrs;
    snr_3layer_list(trial) = -best_val_global;
    
    fprintf('Trial %2d/%2d | SNR_in=%7.3f dB | SNR_1L=%7.3f dB | SNR_3L=%7.3f dB\n', ...
        trial, num_trials, input_snr_list(trial), snr_1layer_list(trial), snr_3layer_list(trial));
end

%% 3. 统计与检测能力评估 ==================================================
mean_snr_in = mean(input_snr_list);
mean_snr_1L = mean(snr_1layer_list);
mean_snr_3L = mean(snr_3layer_list);

std_snr_in = std(input_snr_list);
std_snr_1L = std(snr_1layer_list);
std_snr_3L = std(snr_3layer_list);

gain_1L_vs_in = snr_1layer_list - input_snr_list;
gain_3L_vs_in = snr_3layer_list - input_snr_list;
gain_3L_vs_1L = snr_3layer_list - snr_1layer_list;

% 检测成功率（输出SNR高于输入SNR）
rate_detect_1L = mean(gain_1L_vs_in > 0);
rate_detect_3L = mean(gain_3L_vs_in > 0);

fprintf('\n================= 检测能力统计 =================\n');
fprintf('输入 SNR:      %.4f ± %.4f dB\n', mean_snr_in, std_snr_in);
fprintf('单层输出 SNR:  %.4f ± %.4f dB\n', mean_snr_1L, std_snr_1L);
fprintf('三层输出 SNR:  %.4f ± %.4f dB\n', mean_snr_3L, std_snr_3L);
fprintf('平均增益(1L-输入): %.4f dB\n', mean(gain_1L_vs_in));
fprintf('平均增益(3L-输入): %.4f dB\n', mean(gain_3L_vs_in));
fprintf('平均增益(3L-1L):   %.4f dB\n', mean(gain_3L_vs_1L));
fprintf('检测成功率(1L): %.2f%%\n', 100*rate_detect_1L);
fprintf('检测成功率(3L): %.2f%%\n', 100*rate_detect_3L);
fprintf('================================================\n');

%% 4. 可视化 ===============================================================
figure('Color', 'w', 'Position', [100, 100, 1000, 420]);

subplot(1,2,1);
bar_data = [mean_snr_in, mean_snr_1L, mean_snr_3L];
bar_err  = [std_snr_in, std_snr_1L, std_snr_3L];
b = bar(1:3, bar_data, 0.65); hold on;
set(b, 'FaceColor', [0.35, 0.55, 0.85]);
errorbar(1:3, bar_data, bar_err, 'k.', 'LineWidth', 1.4);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Input', '1-Layer', '3-Layer'});
ylabel('SNR (dB)');
title(sprintf('Alpha噪声下 SNR 对比 (\\alpha=%.1f)', alpha_noise));
grid on;

subplot(1,2,2);
gain_all = [gain_1L_vs_in; gain_3L_vs_in; gain_3L_vs_1L];
group_id = [ones(num_trials,1); 2*ones(num_trials,1); 3*ones(num_trials,1)];
boxplot(gain_all, group_id);
set(gca, 'XTick', 1:3, 'XTickLabel', {'1L-Input', '3L-Input', '3L-1L'});
ylabel('SNR Gain (dB)');
title('检测增益分布 (Monte Carlo)');
grid on;

%% 5. 保存结果 =============================================================
script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, '..', 'results', 'alpha_csr_3layer_detection');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

results.config.fs = fs;
results.config.T = T;
results.config.N = N;
results.config.f0 = f0;
results.config.A0 = A0;
results.config.num_layers = num_layers;
results.config.num_trials = num_trials;
results.config.search_agents = search_agents;
results.config.max_iter = max_iter;
results.config.lb = lb;
results.config.ub = ub;

results.noise.alpha = alpha_noise;
results.noise.beta = beta_noise;
results.noise.sigma = sigma_noise;
results.noise.mu = mu_noise;

results.metrics.input_snr = input_snr_list;
results.metrics.snr_1layer = snr_1layer_list;
results.metrics.snr_3layer = snr_3layer_list;
results.metrics.gain_1L_vs_in = gain_1L_vs_in;
results.metrics.gain_3L_vs_in = gain_3L_vs_in;
results.metrics.gain_3L_vs_1L = gain_3L_vs_1L;
results.metrics.detect_rate_1L = rate_detect_1L;
results.metrics.detect_rate_3L = rate_detect_3L;
results.metrics.layer_snr_all = layer_snr_all;

results.summary.mean_snr_in = mean_snr_in;
results.summary.mean_snr_1L = mean_snr_1L;
results.summary.mean_snr_3L = mean_snr_3L;
results.summary.std_snr_in = std_snr_in;
results.summary.std_snr_1L = std_snr_1L;
results.summary.std_snr_3L = std_snr_3L;

results.params.best_1layer = params_1layer;
results.params.best_3layer = params_3layer;
results.signals.clean_sig = clean_sig;
results.signals.input_noisy_all = input_noisy_all;
results.signals.layer_output_all = layer_output_all;

save(fullfile(result_dir, 'res_alpha_csr_3layer_detection.mat'), 'results');

fprintf('结果已保存到: %s\n', fullfile(result_dir, 'res_alpha_csr_3layer_detection.mat'));


%% ========================= 辅助函数 ======================================
function fitness = CostFunc_HSUBSR_Layer(theta, clean_sig, noise_force, fs, f0)
% 单层 HSUBSR 在给定输入下的适应度（最小化 -SNR）
try
    [a, b, k1, k2] = CalibrateHSUBSR(theta(1), theta(2), theta(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    x_out = RK4Solver2(drift, clean_sig, noise_force, fs);
    x_out = x_out - mean(x_out);
    
    idx0 = round(0.1 * length(x_out));
    x_steady = x_out(idx0+1:end);
    
    snr_val = SNRo2(x_steady, fs, f0);
    if isnan(snr_val) || isinf(snr_val)
        fitness = 1e5;
        return;
    end
    
    fitness = -snr_val;
catch
    fitness = 1e5;
end
end

function fitness = CostFunc_HSUBSR_Global(theta_vec, clean_sig, noise_force, fs, f0, num_layers)
% 三层级联 HSUBSR 全局参数适应度（最小化 -SNR）
try
    current_clean = clean_sig;
    current_noise = noise_force;
    n_samples = length(clean_sig);
    
    for k = 1:num_layers
        idx_start = (k-1)*3 + 1;
        theta_k = theta_vec(idx_start:idx_start+2);
        
        [a, b, k1, k2] = CalibrateHSUBSR(theta_k(1), theta_k(2), theta_k(3));
        drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
        
        x_out = RK4Solver2(drift, current_clean, current_noise, fs);
        x_out = x_out - mean(x_out);
        
        current_clean = x_out;
        current_noise = zeros(n_samples, 1); % 后续层不再注入新噪声
    end
    
    idx0 = round(0.1 * length(current_clean));
    x_steady = current_clean(idx0+1:end);
    
    snr_val = SNRo2(x_steady, fs, f0);
    if isnan(snr_val) || isinf(snr_val)
        fitness = 1e5;
        return;
    end
    
    fitness = -snr_val;
catch
    fitness = 1e5;
end
end

function [layer_outputs, layer_snrs] = SimulateHSUBSRCascade(theta_vec, clean_sig, noise_force, fs, f0, num_layers)
% 在给定全局参数下仿真级联系统，返回每层输出及每层SNR
n_samples = length(clean_sig);
layer_outputs = zeros(n_samples, num_layers);
layer_snrs = zeros(1, num_layers);

current_clean = clean_sig;
current_noise = noise_force;

for k = 1:num_layers
    idx_start = (k-1)*3 + 1;
    theta_k = theta_vec(idx_start:idx_start+2);
    
    [a, b, k1, k2] = CalibrateHSUBSR(theta_k(1), theta_k(2), theta_k(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    x_out = RK4Solver2(drift, current_clean, current_noise, fs);
    x_out = x_out - mean(x_out);
    
    layer_outputs(:, k) = x_out;
    
    idx0 = round(0.1 * n_samples);
    x_steady = x_out(idx0+1:end);
    layer_snrs(k) = SNRo2(x_steady, fs, f0);
    
    current_clean = x_out;
    current_noise = zeros(n_samples, 1);
end
end

function x = GenerateAlphaStableNoise(alpha, beta, sigma, mu, m, n)
% CMS 算法生成 Alpha 稳定分布噪声
if nargin < 6, n = 1; end
if nargin < 5, m = 1; end
if nargin < 4, mu = 0; end
if nargin < 3, sigma = 1; end

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
