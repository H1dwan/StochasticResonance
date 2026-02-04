% =========================================================================
% Description: 对比固定参数寻优、逐层参数寻优和全局参数寻优，单次实验
%
% Author: LiuShuang
% Created: 2026-02-03
% Last Modified: 2026-02-03
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;
rng(0);  % 固定随机数种子，便于结果复现

%% 1. 仿真信号与环境设置 ====================================================
fs = 5;              % 采样频率 Hz
T  = 2000;           % 信号长度
N  = fs * T;
t  = (0:N-1)'/fs;

% 目标信号
f0 = 0.01;            % 目标频率
A0 = 0.10;            % 信号幅值
clean_sig = A0 * sin(2*pi*f0*t);

% 背景噪声
D  = 0.4;             % 噪声强度
h  = 1/fs;
noise_seq = sqrt(2*D/h) * randn(N, 1);

input_sig = clean_sig + noise_seq;
snr_in = SNRo2(input_sig, fs, f0);

fprintf('>>> 仿真环境配置:\n');
fprintf('    采样率: %d Hz, 时长: %d s\n', fs, T);
fprintf('    输入 SNR: %.4f dB\n', snr_in);


%% 2. PSO 参数设置 ====================================================
num_layers = 3;      % 级联层数
search_agents   = 20; % PSO 种群规模
max_iter        = 50; % PSO 迭代次数
lb = [0.5, 0.15, 1.001];  % 单组参数下界 [xm, dU, shape]
ub = [2.0, 1.0, 1.9]; % 单组参数上界 [xm, dU, shape]

%% 3. 三种策略对比 ====================================================

% 策略 A: 固定参数寻优 (Fixed Parameter Optimization)
% 假设所有层使用相同的 [xm, dU, shape]
fobj_fixed = @(theta) CostFunc_Fixed(theta, input_sig, fs, f0, num_layers);
[best_val_fixed, best_pos_fixed, ~] = PSO(search_agents, max_iter, lb, ub, length(lb), fobj_fixed);
snr_fixed = -best_val_fixed;    % 适应度取负还原为 SNR
params_fixed = best_pos_fixed;

% 策略 B: 逐层寻优 (Layer-by-Layer Optimization)
current_input = input_sig;
current_noise = zeros(size(input_sig)); % 初始噪声包含在 input_sig 中，后续层视为无噪输入

layer_params = zeros(1, 3 * num_layers);
for lvl = 1:num_layers
    % 定义单层适应度函数 (利用现有的 SNREvaluator)
    % 注意: SNREvaluator 内部使用了 CalibrateHSUBSR，符合 HSUBSR 模型要求
    fobj_layer = @(theta) SNREvaluator(theta, current_input, current_noise, fs, f0);
    
    % 执行 PSO
    [~, best_pos_layer, ~] = PSO(search_agents, max_iter, lb, ub, 3, fobj_layer);
    
    % 记录该层最优参数
    layer_params((lvl-1)*3+1:(lvl-1)*3+3) = best_pos_layer;
    
    % 使用最优参数生成该层输出，作为下一层输入
    % 参数映射: [xm, dU, shape] -> [a, b, k1, k2]
    [a, b, k1, k2] = CalibrateHSUBSR(best_pos_layer(1), best_pos_layer(2), best_pos_layer(3));
    drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    current_input = RK4Solver(drift_func, current_input + current_noise, 1/fs);
    
    % 去直流
    current_input = current_input - mean(current_input);
    
    % 第一层处理完后，后续层的输入本身已包含噪声，不再叠加额外噪声
    current_noise = zeros(size(current_input));
end

% 计算最终层 SNR
% 截取稳态数据 (后90%)
steady_len = round(0.9 * length(current_input));
sig_steady = current_input(end-steady_len+1:end);
snr_layer = SNRo2(sig_steady, fs, f0);

params_layer = layer_params;

% 策略 C: 全局寻优 (Global Optimization)
dim_global = 3 * num_layers;
lb_global = repmat(lb, 1, num_layers); % 扩展下界
ub_global = repmat(ub, 1, num_layers); % 扩展上界

fobj_global = @(theta) CostFunc_Global(theta, input_sig, fs, f0, num_layers);

[best_val_global, best_pos_global, ~] = PSO(search_agents, max_iter, lb_global, ub_global, dim_global, fobj_global);
snr_global = -best_val_global;
params_global = best_pos_global;


%% 辅助函数 ====================================================
function fitness = CostFunc_Fixed(theta_vec, input_sig, fs, f0, levels)
% CostFunc_Fixed: 计算固定参数下的级联系统适应度
% theta_vec: [xm, dU, shape] (仅一组参数，应用于所有层)

% 1. 参数校准
[a, b, k1, k2] = CalibrateHSUBSR(theta_vec(1), theta_vec(2), theta_vec(3));
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);

current_in = input_sig;

% 2. 级联循环
for i = 1:levels
    % 每一层都使用相同的 drift_func
    % RK4Solver 求解: dx = (-U'(x) + input) * dt
    current_in = RK4Solver(drift_func, current_in, 1/fs);
    
    % 级联处理关键: 去除直流分量，防止信号漂移出势阱
    current_in = current_in - mean(current_in);
end

% 3. 计算末级 SNR
steady_idx = round(0.1 * length(current_in));
sig_steady = current_in(steady_idx+1:end);

val = SNRo2(sig_steady, fs, f0);
if isinf(val), val = -100; end % 异常处理
fitness = -val; % 最小化目标
end

function fitness = CostFunc_Global(theta_vec, input_sig, fs, f0, levels)
% CostFunc_Global: 计算全局参数下的级联系统适应度
% theta_vec: [xm1, dU1, shape1, xm2, dU2, shape2, ...] (长度 = 3 * levels)

current_in = input_sig;

% 1. 级联循环
for i = 1:levels
    % 提取当前层的参数索引
    idx_start = (i-1)*3 + 1;
    p = theta_vec(idx_start : idx_start+2);
    
    % 参数校准
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    % 求解
    current_in = RK4Solver(drift_func, current_in, 1/fs);
    current_in = current_in - mean(current_in);
end

% 2. 计算末级 SNR
steady_idx = round(0.1 * length(current_in));
sig_steady = current_in(steady_idx+1:end);

val = SNRo2(sig_steady, fs, f0);
if isinf(val), val = -100; end
fitness = -val;
end