% =========================================================================
% Description: 验证 CSR 在 CWRU 轴承数据上的性能
%
% Author: LiuShuang
% Created: 2026-02-10
% Last Modified: 2026-02-10
%
% Usage: 保证路径下存在 CWRU 数据文件
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;
% seed = 28;
% rng(seed);

%% 1. 加载数据 ==============================================================
data_filename = '135.mat';
S = load(data_filename);  % 加载到结构体

% 查找变量名
var_names = fieldnames(S);
de_idx = find(contains(var_names, '_DE_time'));
rpm_idx = find(contains(var_names, 'RPM'));

if isempty(de_idx) || isempty(rpm_idx)
    error('找不到所需的变量');
end

% 提取数据
raw_sig = S.(var_names{de_idx(1)});
rpm = S.(var_names{rpm_idx(1)});
fs_raw = 48000;      % CWRU标准采样率 12kHz

% 计算轴承故障频率
[freqs, desc] = BearingHz(rpm);
f_fault_real = freqs.BPFO;      % 内/外圈故障频率
% f_fault_real = freqs.BPFI;      % 内/外圈故障频率

%% 2. 信号预处理 ============================================================
% 截取数据片段进行分析（避免计算量过大）
N_sample = 16384;
N = length(raw_sig);
t = (0:N-1)/fs_raw;
sig_segment = raw_sig(1:N_sample);

% 输入信号的频谱分析
fprintf('Input snr: %.4fdb\n', SNRo2(sig_segment, fs_raw, f_fault_real));
Plot_Time_Frequency2(sig_segment, fs_raw, 'FaultFreq', f_fault_real);

% 包络提取，SR作为低通滤波器，直接处理高频载波效果差，需先提取故障包络
sig_env = abs(hilbert(sig_segment));

% 2.2 去直流
sig_ac = sig_env - mean(sig_env);

% 2.3 幅值尺度变换 (Amplitude Rescaling)
% z-score 标准化
sig_mean = mean(sig_ac);
sig_std = std(sig_ac);
sig_processed = (sig_ac - sig_mean) / sig_std;

% 2.4 频率尺度变换 (Frequency Rescaling)
% 将真实采样频率映射到固定的5Hz
R = 9600;
fs = fs_raw / R;    % 变换后的采样频率
f_target_sr = f_fault_real / R;     % 变换后的目标频率

% 更新数值积分步长 h
% 原始步长 h0 = 1/fs. 变换后系统"感知"到的步长需放大 R 倍
% 理论依据: d x / d (t/R) = ... => 相当于时间变慢，频率变低
h_sr = (1/fs_raw) * R;

fprintf('frequency rescaling factor R: %.2f\n', R);
fprintf('Original target frequency: %.3f Hz, Rescaled target frequency: %.3f Hz\n', f_fault_real, f_target_sr);
fprintf('Original fs: %.2f Hz, Rescaled fs: %.2f Hz\n', fs_raw, fs);

%% 3. 自适应 PSO 寻优 ==============================

% PSO 参数设置
dim = 3;              % 优化变量维度 [xm, dU, shape]
lb = [0.5, 0.1, 1.01];  % 单组参数下界 [xm, dU, shape]
ub = [3.0, 2.0, 1.99]; % 单组参数上界 [xm, dU, shape]

SearchAgents_no = 20; % 种群规模
Max_iter = 50;        % 迭代次数

% 多次全局优化，取最优解以降低随机性
n_pso_runs = 1;

num_layers = 3;      % 级联层数

% 全局寻优 (Global Optimization)
dim_global = 3 * num_layers;
lb_global = repmat(lb, 1, num_layers); % 扩展下界
ub_global = repmat(ub, 1, num_layers); % 扩展上界

fobj_global = @(theta) CostFunc_Global(theta, sig_processed, fs, f_target_sr, num_layers);

fprintf('\nStarting PSO Optimization...\n');
best_val_global = inf;
best_pos_global = [];
best_curve_global = [];

for run_id = 1:n_pso_runs
    [val_tmp, pos_tmp, curve_tmp] = PSO(SearchAgents_no, Max_iter, lb_global, ub_global, dim_global, fobj_global, false);
    if val_tmp < best_val_global
        best_val_global = val_tmp;
        best_pos_global = pos_tmp;
        best_curve_global = curve_tmp;
    end
    fprintf('  Run %d/%d -> Best SNR: %.4f dB\n', run_id, n_pso_runs, -val_tmp);
end

snr_global = -best_val_global;
params_global = best_pos_global;

fprintf('Global optimization finished.\n');
fprintf('Best 3-layer global SNR (SR domain): %.4f dB\n', snr_global);
fprintf('Best params: %s\n', mat2str(params_global, 4));


%% 4. 最优结果验证与可视化 ============================================================

% 4.1 用全局最优参数重建三级级联输出
[layer_outputs, layer_snr_sr] = SimulateCSR_Global(params_global, sig_processed, fs, f_target_sr, num_layers);
x_final = layer_outputs(:, end);

% 4.2 输入与各层性能统计
snr_input_sr = SNRo2(sig_processed, fs, f_target_sr);

fprintf('\n=========== CWRU: 三级级联全局寻优性能 ===========\n');
fprintf('Input SNR (SR domain): %.4f dB\n', snr_input_sr);
for k = 1:num_layers
    fprintf('Layer-%d SNR (SR domain): %.4f dB\n', k, layer_snr_sr(k));
end
fprintf('Final Gain vs Input: %.4f dB\n', layer_snr_sr(end) - snr_input_sr);
fprintf('===============================================\n');

% 4.3 可视化（时域）
t_sr = (0:length(sig_processed)-1) / fs;
show_len = length(sig_processed);
idx_show = 1:show_len;

figure('Color', 'w', 'Position', [100, 80, 1000, 780]);
subplot(num_layers+1, 1, 1);
plot(t_sr(idx_show), sig_processed(idx_show), 'k', 'LineWidth', 1.0);
grid on;
ylabel('Amp');
title(sprintf('Input (SR domain), SNR=%.3f dB', snr_input_sr));

for k = 1:num_layers
    subplot(num_layers+1, 1, k+1);
    plot(t_sr(idx_show), layer_outputs(idx_show, k), 'LineWidth', 1.0);
    grid on;
    ylabel('Amp');
    title(sprintf('Layer-%d Output, SNR=%.3f dB', k, layer_snr_sr(k)));
end
xlabel('Time (s)');

% 4.4 可视化（最终输出频域，映射回原频率）
fs_restore = fs * R;
f_fault_restore = f_target_sr * R;
x_final = x_final - mean(x_final); % 去直流分量
x_final = 1 * x_final(0.1*end:end);
Plot_Time_Frequency2(x_final, fs_restore, 'FaultFreq', f_fault_restore, 'FaultLabel', '$f_{BPFO}=$');

% 4.5 保存结果
script_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(script_dir, '..', 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

[~, data_tag, ~] = fileparts(data_filename);
save_name = sprintf('res_%s_2026_csr_global.mat', data_tag);

results.name = 'HSUBSR-CSR-Global-3Layer';
results.config.data_filename = data_filename;
% results.config.seed = seed;
results.config.fs_raw = fs_raw;
results.config.fs_sr = fs;
results.config.rescale_R = R;
results.config.f_fault_real = f_fault_real;
results.config.f_target_sr = f_target_sr;
results.config.num_layers = num_layers;
results.config.search_agents = SearchAgents_no;
results.config.max_iter = Max_iter;
results.config.n_pso_runs = n_pso_runs;
results.config.lb = lb;
results.config.ub = ub;

results.best.best_val = best_val_global;
results.best.best_snr_sr = snr_global;
results.best.best_params = params_global;
results.best.best_curve = best_curve_global;

results.metrics.input_snr_sr = snr_input_sr;
results.metrics.layer_snr_sr = layer_snr_sr;
results.metrics.final_gain_vs_input = layer_snr_sr(end) - snr_input_sr;

results.signals.input_sr = sig_processed;
results.signals.layer_outputs_sr = layer_outputs;
results.signals.final_output_sr = x_final;

save(fullfile(result_dir, save_name), 'results');
fprintf('Results saved: %s\n', fullfile(result_dir, save_name));


%% 辅助函数
function [freqs, description] = BearingHz(Fr_rpm, D, d, Z, alpha_deg, verbose)
% BearingHz 计算滚动轴承的故障特征频率
%
% Useage:
%   freqs = BearingHz(1797); % 使用 CWRU 默认参数 (SKF 6205-2RS)
%   freqs = BearingHz(1797, 39.04, 7.94, 9, 0); % 自定义参数
%
% Input:
%   Fr_rpm    : 轴旋转速度 (RPM), 支持向量输入
%   D         : 节径 (Pitch Diameter), 单位 mm (默认: CWRU参数 39.04mm)
%   d         : 滚动体直径 (Ball Diameter), 单位 mm (默认: CWRU参数 7.94mm)
%   Z         : 滚动体个数 (Number of Elements) (默认: 9)
%   alpha_deg : 接触角 (Contact Angle), 单位 度 (默认: 0)
%   verbose   : 是否打印结果 (true/false, 默认 true)
%
% Output:
%   freqs     : 包含各故障频率的结构体 (Hz)
%       .BPFI : 内圈故障频率
%       .BPFO : 外圈故障频率
%       .BSF  : 滚动体故障频率
%       .FTF  : 保持架(基频)频率
%       .Fr   : 转频
%
% Reference:
%   CWRU Bearing: SKF 6205-2RS JEM SKF
%   D = 1.537 inch ≈ 39.04 mm
%   d = 0.3126 inch ≈ 7.94 mm
%   Z = 9 balls
%   alpha = 0 deg

% =========================================================================
% 1. 参数预处理与默认值设置 (针对 CWRU 数据集)
% =========================================================================
if nargin < 6, verbose = true; end
% 如果只输入了转速，默认加载 SKF 6205-2RS 参数
if nargin < 5 || (isempty(D) && isempty(d))
    % SKF 6205-2RS 参数 (CWRU标准)
    D = 39.04;
    d = 7.94;
    Z = 9;
    alpha_deg = 0;
    if verbose && nargin < 2
        fprintf('Info: Using default CWRU (SKF 6205-2RS) parameters.\n');
    end
end

% 角度转弧度
alpha = alpha_deg * pi / 180;

% 转速 RPM -> Hz
fr = Fr_rpm / 60;

% 几何比率 (简化公式书写)
ratio = (d / D) * cos(alpha);

% =========================================================================
% 2. 频率计算 (Hz)
% =========================================================================

% 内圈故障频率 (Inner Race)
bpfi = (Z / 2) * (1 + ratio) .* fr;

% 外圈故障频率 (Outer Race)
bpfo = (Z / 2) * (1 - ratio) .* fr;

% 保持架故障频率 (Cage / Fundamental Train)
ftf  = (1 / 2) * (1 - ratio) .* fr;

% 滚动体故障频率 (Ball Spin)
% 公式: D/(2d) * (1 - (d/D*cos(a))^2) * fr
bsf  = (D / (2 * d)) * (1 - ratio.^2) .* fr;

% =========================================================================
% 3. 结果封装与显示
% =========================================================================
freqs.Fr   = fr;
freqs.BPFI = bpfi;
freqs.BPFO = bpfo;
freqs.BSF  = bsf;
freqs.FTF  = ftf;

% 生成描述文本（可选）
description = sprintf(['Rotation Frequencies (Fr = %.2f Hz):\n' ...
    '  Inner Race (BPFI): %.2f Hz (%.2f X)\n' ...
    '  Outer Race (BPFO): %.2f Hz (%.2f X)\n' ...
    '  Ball Spin  (BSF) : %.2f Hz (%.2f X)\n' ...
    '  Cage       (FTF) : %.2f Hz (%.2f X)'], ...
    mean(fr), ...
    mean(bpfi), mean(bpfi)/mean(fr), ...
    mean(bpfo), mean(bpfo)/mean(fr), ...
    mean(bsf),  mean(bsf)/mean(fr), ...
    mean(ftf),  mean(ftf)/mean(fr));

if verbose
    fprintf('--------------------------------------------------\n');
    fprintf('%s\n', description);
    fprintf('--------------------------------------------------\n');
end

end

% --- 1. UBSR 映射与适应度 ---
function [a, b] = MapParams_UBSR(xm, dU)
% 物理推导: U(x) = -a/2 x^2 + b/4 x^4
% xm^2 = a/b; dU = a^2/(4b)
% 解得: a = 4*dU / xm^2; b = 4*dU / xm^4
a = 4 * dU / (xm^2);
b = 4 * dU / (xm^4);
end

function val = Fitness_UBSR_Fair(p, s, h, f0)
[a, b] = MapParams_UBSR(p(1), p(2)); % p=[xm, dU]
drift = @(x) UBSR_Dynamics(x, a, b);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end):end);
val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_UBSR_Fair(p, s, h)
[a, b] = MapParams_UBSR(p(1), p(2));
drift = @(x) UBSR_Dynamics(x, a, b);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end)+1:end);
end

% --- 2. PLBSR 映射与适应度 ---
function [U0, L0] = MapParams_PLBSR(xm, dU)
% PLBSR 参数直接对应: L0是阱宽(xm), U0是垒高(dU)
L0 = xm;
U0 = dU;
end

function val = Fitness_PLBSR_Fair(p, s, h, f0)
[U0, L0] = MapParams_PLBSR(p(1), p(2)); % p=[xm, dU]
drift = @(x) PLBSR_Dynamics(x, U0, L0);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end):end);
val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_PLBSR_Fair(p, s, h)
[U0, L0] = MapParams_PLBSR(p(1), p(2));
drift = @(x) PLBSR_Dynamics(x, U0, L0);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end)+1:end);
end

% --- 3. HSUBSR 映射与适应度 ---
function val = Fitness_HSUBSR_Fair(p, s, h, f0)
% p = [xm, dU, shape]
% 调用已有的校准函数
[a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end):end);
val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_HSUBSR_Fair(p, s, h)
[a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x = RK4Solver(drift, s, h);
x = x(round(0.1*end)+1:end);
end

function Plot_Time_Frequency(x, fs, options)
% 绘制实信号的单边幅频图
%   该函数用于同时绘制信号的时域波形和频域幅频特性曲线
%
% 输入参数:
%   x - 输入实信号，一维数组
%   fs - 采样频率，单位Hz
%   N - 信号长度，即信号点数
%   options - 绘图选项结构体
%       .LineWidth - 线条宽度，默认为1

% 参数解析
arguments
    x   double
    fs  double
    options.LineWidth (1,1) {mustBeNumeric} = 1
end

N = length(x);
% 时间轴和频率轴计算
t = (0:N-1) / fs;
f = fs/N*(0:(N/2)); % 频率范围（包含奈奎斯特频率）

% FFT计算和幅值归一化
P2 = abs(fft(x)/N); % 计算fft并归一化纵轴（双边）
P1 = P2(1:N/2+1);   % 截取单边谱
P1(2:end-1) = 2*P1(2:end-1);    % 能量还原（除0Hz和奈奎斯特频率）

% 绘制时域和频域图形
SetThesisDefaultStyle();
CreateThesisFigure(8, 3);
layout = tiledlayout(1,2);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
xlabel('Time[s]')
ylabel('Amplitude')
% xlim([0 0.3])

nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlim([0 1000])
xlabel('Frequency[Hz]')
ylabel('Amplitude')

% Mark the highest spectral peak (exclude DC)
ax = gca;
if numel(P1) > 1
    xlimv = ax.XLim;
    ylimv = ax.YLim;
    inRange = f >= xlimv(1) & f <= xlimv(2);
    if any(inRange)
        [peakVal, relIdx] = max(P1(inRange));
        idxs = find(inRange);
        peakIdx = idxs(relIdx);
        peakFreq = f(peakIdx);
        axpos = ax.Position; % normalized in figure
        xnorm = axpos(1) + (peakFreq - xlimv(1)) / diff(xlimv) * axpos(3);
        ynorm = axpos(2) + (peakVal - ylimv(1)) / diff(ylimv) * axpos(4) + 0.01;
        % Clamp to keep annotation within figure bounds
        x1 = min(max(xnorm, 0.02), 0.98);
        y1 = min(max(ynorm, 0.02), 0.98);
        x0 = min(max(x1 + 0.06, 0.02), 0.98);
        y0 = min(max(y1 - 0.06, 0.02), 0.98);
        annotation('textarrow', [x0 x1], [y0 y1], ...
            'String', sprintf('$f=%.2f\\,\\mathrm{Hz}$', peakFreq), ...
            'FontSize', 12, 'LineWidth', 1, 'Interpreter', 'latex');
    end
end

end

function Plot_Time_Frequency2(x, fs, options)
% 绘制实信号的单边幅频图
%   该函数用于同时绘制信号的时域波形和频域幅频特性曲线
%
% 输入参数:
%   x - 输入实信号，一维数组
%   fs - 采样频率，单位Hz
%   N - 信号长度，即信号点数
%   options - 绘图选项结构体
%       .LineWidth - 线条宽度，默认为1

% 参数解析
arguments
    x   double
    fs  double
    options.LineWidth (1,1) {mustBeNumeric} = 1
    options.FaultFreq (1,1) double = NaN
    options.FaultLabel (1,1) string = '$f_{BPFO}=$'
end

N = length(x);
% 时间轴和频率轴计算
t = (0:N-1) / fs;
f = fs/N*(0:(N/2)); % 频率范围（包含奈奎斯特频率）

% FFT计算和幅值归一化
P2 = abs(fft(x)/N); % 计算fft并归一化纵轴（双边）
P1 = P2(1:N/2+1);   % 截取单边谱
P1(2:end-1) = 2*P1(2:end-1);    % 能量还原（除0Hz和奈奎斯特频率）

% 绘制时域和频域图形
SetThesisDefaultStyle();
CreateThesisFigure(8, 6);
layout = tiledlayout(2,1);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
xlabel('Time[s]')
ylabel('Amplitude')
xlim([0 0.3])

ax = nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlabel('Frequency[Hz]')
ylabel('Amplitude')
xlim([0 1000])


if ~isnan(options.FaultFreq)
    % ff = options.FaultFreq;
    [~, ff_idx] = max(P1);
    ff = f(ff_idx); % 直接标记频谱峰值对应的频率
    yval = interp1(f, P1, ff, 'linear', 'extrap');
    yval = max(0, min(yval, max(P1)));
    
    axpos = ax.Position; % normalized in figure
    xnorm = axpos(1) + (ff - ax.XLim(1)) / diff(ax.XLim) * axpos(3);
    ynorm = axpos(2) + (yval - ax.YLim(1)) / diff(ax.YLim) * axpos(4);
    
    % 目标点(箭头尖端)与起点均约束到[0,1]，并固定方向为右下->左上(↖)
    x_tip = min(max(xnorm, 0.02), 0.98);
    y_tip = min(max(ynorm, 0.02), 0.98);
    
    x_tail = min(max(x_tip + 0.06, 0.02), 0.98); % 更靠右
    y_tail = min(max(y_tip - 0.06, 0.02), 0.98); % 更靠下
    
    annotation('textarrow', [x_tail x_tip], [y_tail y_tip], ...
        'String', sprintf('%s %.2f Hz', options.FaultLabel, ff), ...
        'FontSize', 12, 'LineWidth', 1, 'Interpreter', 'latex');
end

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

function [layer_outputs, layer_snrs] = SimulateCSR_Global(theta_vec, input_sig, fs, f0, levels)
% SimulateCSR_Global: 在给定全局参数下重建各层输出与层级SNR
n_samples = length(input_sig);
layer_outputs = zeros(n_samples, levels);
layer_snrs = zeros(1, levels);

current_in = input_sig;
for i = 1:levels
    idx_start = (i-1)*3 + 1;
    p = theta_vec(idx_start:idx_start+2);
    
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    current_in = RK4Solver(drift_func, current_in, 1/fs);
    current_in = current_in - mean(current_in);
    
    layer_outputs(:, i) = current_in;
    
    steady_idx = round(0.1 * n_samples);
    sig_steady = current_in(steady_idx+1:end);
    layer_snrs(i) = SNRo2(sig_steady, fs, f0);
end
end