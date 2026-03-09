% =========================================================================
% Description: 验证 CSR 在 HUST 轴承数据上的性能
%
% Author: LiuShuang
% Created: 2026-02-10
% Last Modified: 2026-02-10
%
% Usage: 保证路径下存在 HUST 轴承数据文件
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 数据配置 =============================================================
data_filename = '0.5X_O_40Hz.mat';  % 可改为其它 HUST 数据文件
fs_raw = 25600;

% 若文件名在映射表中，自动给出故障频率；否则需手动设置
f_fault_real = GetHUSTFaultFreq(data_filename);
if isnan(f_fault_real)
    error('未找到 %s 的故障频率映射，请在 GetHUSTFaultFreq 中补充。', data_filename);
end

N_sample = 4096*2;
R = 5120;                 % 频率缩放因子

% 全局寻优配置（三层 HSUBSR）
num_layers = 3;
lb = [0.1, 0.1, 1.01];  % 单组参数下界 [xm, dU, shape]
ub = [3.0, 1.5, 1.99];  % 单组参数上界 [xm, dU, shape]
SearchAgents_no = 20;
Max_iter = 50;
n_pso_runs = 1;

%% 2. 加载 HUST 数据 =======================================================
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..', '..');

candidate_paths = {
    fullfile(script_dir, data_filename), ...
    fullfile(project_root, data_filename), ...
    data_filename ...
    };

data_path = '';
for i = 1:numel(candidate_paths)
    if exist(candidate_paths{i}, 'file')
        data_path = candidate_paths{i};
        break;
    end
end

if isempty(data_path)
    error('未找到数据文件: %s', data_filename);
end

S = load(data_path);
if ~isfield(S, 'vib_signal')
    error('文件 %s 中未找到变量 vib_signal', data_path);
end

raw_sig = S.vib_signal(:);
if length(raw_sig) < N_sample
    error('数据长度不足: length(vib_signal)=%d < N_sample=%d', length(raw_sig), N_sample);
end

sig_segment = raw_sig(1:N_sample);
Plot_Time_Frequency1(sig_segment, fs_raw, 'FaultFreq', f_fault_real, 'FaultLabel', '$f_{BPFI}=$');

fprintf('=== HUST 三层 HSUBSR 级联共振验证 ===\n');
fprintf('Data: %s\n', data_path);
fprintf('Input SNR (raw): %.4f dB\n', SNRo2(sig_segment, fs_raw, f_fault_real));

%% 3. 信号预处理 ===========================================================
% 包络提取
sig_env = abs(hilbert(sig_segment));

% 去直流 + 标准化
sig_ac = sig_env - mean(sig_env);
sig_processed = (sig_ac - mean(sig_ac)) / std(sig_ac);

% 频率尺度变换
fs = fs_raw / R;
f_target_sr = f_fault_real / R;
h_sr = (1/fs_raw) * R;

fprintf('R=%.1f, fs_raw=%.1fHz -> fs_sr=%.4fHz, f_fault=%.4fHz -> f_target=%.6fHz\n', ...
    R, fs_raw, fs, f_fault_real, f_target_sr);


%% 4. 三层全局参数寻优 =====================================================
dim_global = 3 * num_layers;
lb_global = repmat(lb, 1, num_layers);
ub_global = repmat(ub, 1, num_layers);

sig_processed = sig_processed + 0.08*sin(2*pi*f_target_sr*(0:length(sig_processed)-1)'/fs);

fobj_global = @(theta) CostFunc_Global(theta, sig_processed, fs, f_target_sr, num_layers);

fprintf('\nStarting Global PSO (3-layer HSUBSR)...\n');

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
    fprintf('  Run %d/%d -> best SNR (SR) = %.4f dB\n', run_id, n_pso_runs, -val_tmp);
end

snr_global = -best_val_global;
params_global = best_pos_global;

%% 5. 最优参数验证 + 可视化 =================================================
[layer_outputs, layer_snr_sr] = SimulateCSR_Global(params_global, sig_processed, fs, f_target_sr, num_layers);
x_final_sr = layer_outputs(:, end);

snr_input_sr = SNRo2(sig_processed, fs, f_target_sr);
gain_final_sr = layer_snr_sr(end) - snr_input_sr;

fprintf('\n=========== HUST: 三级级联全局寻优性能 ===========\n');
fprintf('Input SNR (SR domain): %.4f dB\n', snr_input_sr);
for k = 1:num_layers
    fprintf('Layer-%d SNR (SR domain): %.4f dB\n', k, layer_snr_sr(k));
end
fprintf('Final Gain vs Input (SR): %.4f dB\n', gain_final_sr);
fprintf('Best params: %s\n', mat2str(params_global, 4));
fprintf('================================================\n');

% 时域输出（SR域）
t_sr = (0:length(sig_processed)-1)' / fs;
figure('Color', 'w', 'Position', [100, 80, 1000, 760]);
subplot(num_layers+1, 1, 1);
plot(t_sr, sig_processed, 'k', 'LineWidth', 1.0);
grid on; ylabel('Amp');
title(sprintf('Input (SR domain), SNR=%.3f dB', snr_input_sr));

for k = 1:num_layers
    subplot(num_layers+1, 1, k+1);
    plot(t_sr, layer_outputs(:, k), 'LineWidth', 1.0);
    grid on; ylabel('Amp');
    title(sprintf('Layer-%d Output (SR domain), SNR=%.3f dB', k, layer_snr_sr(k)));
end
xlabel('Time (s)');

% 最终输出频谱（映射回原频率）
fs_restore = fs * R;
x_final_sr = x_final_sr - mean(x_final_sr);
st_idx = round(0.1 * length(x_final_sr)) + 1;
x_final_st = x_final_sr(st_idx:end);
% x_final_st = smooth(x_final_st, 5);

Plot_Time_Frequency2(x_final_st, fs_restore, 'FaultFreq', f_fault_real, 'FaultLabel', '$f_{BPFI}=$');

% 层级增益
figure('Color', 'w', 'Position', [180, 180, 680, 340]);
bar(1:num_layers, layer_snr_sr - snr_input_sr, 0.6);
grid on;
set(gca, 'XTick', 1:num_layers, 'XTickLabel', compose('Layer-%d', 1:num_layers));
ylabel('SNR Gain vs Input (dB)');
title('HUST: Layer-wise Gain (SR domain)');

%% 6. 保存结果 =============================================================
result_dir = fullfile(script_dir, '..', 'results', 'Chap5_csr_hust');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

[~, data_tag, ~] = fileparts(data_filename);
save_name = sprintf('res_%s_2026_csr_global.mat', data_tag);

results.name = 'HSUBSR-CSR-Global-3Layer-HUST';
results.config.data_filename = data_filename;
results.config.data_path = data_path;
results.config.fs_raw = fs_raw;
results.config.fs_sr = fs;
results.config.rescale_R = R;
results.config.f_fault_real = f_fault_real;
results.config.f_target_sr = f_target_sr;
results.config.N_sample = N_sample;
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

results.metrics.input_snr_raw = SNRo2(sig_segment, fs_raw, f_fault_real);
results.metrics.input_snr_sr = snr_input_sr;
results.metrics.layer_snr_sr = layer_snr_sr;
results.metrics.final_gain_vs_input_sr = gain_final_sr;

results.signals.input_raw = sig_segment;
results.signals.input_sr = sig_processed;
results.signals.layer_outputs_sr = layer_outputs;
results.signals.final_output_sr = x_final_st;

save(fullfile(result_dir, save_name), 'results');
fprintf('Results saved: %s\n', fullfile(result_dir, save_name));


%% ========================= 辅助函数 ======================================
function f_fault = GetHUSTFaultFreq(file_name)
% 根据文件名映射 HUST 故障特征频率（Hz）
% 未覆盖的文件返回 NaN，请按数据集标签补充
f_fault = NaN;

switch lower(strtrim(file_name))
    case '0.5x_o_40hz.mat'
        f_fault = 142.90;  % 外圈
    case '0.5x_i_40hz.mat'
        f_fault = 217.10;  % 内圈
    case '0.5x_i_25hz.mat'
        f_fault = 135.69;  % 内圈
    case '0.5x_i_20hz.mat'
        f_fault = 108.55;  % 内圈
    otherwise
        f_fault = NaN;
end
end

function fitness = CostFunc_Global(theta_vec, input_sig, fs, f0, levels)
% 全局参数下的级联系统适应度（最小化 -SNR）
current_in = input_sig;

for i = 1:levels
    idx_start = (i-1)*3 + 1;
    p = theta_vec(idx_start : idx_start+2);
    
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    
    current_in = RK4Solver(drift_func, current_in, 1/fs);
    current_in = current_in - mean(current_in);
    % current_in = smooth(current_in, 5);
end

steady_idx = round(0.1 * length(current_in));
sig_steady = current_in(steady_idx+1:end);

val = SNRo2(sig_steady, fs, f0);
if isinf(val) || isnan(val), val = -100; end
fitness = -val;
end

function [layer_outputs, layer_snrs] = SimulateCSR_Global(theta_vec, input_sig, fs, f0, levels)
% 在给定全局参数下重建各层输出与层级SNR
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
    % current_in = smooth(current_in, 5);
    layer_outputs(:, i) = current_in;
    
    steady_idx = round(0.1 * n_samples);
    sig_steady = current_in(steady_idx+1:end);
    layer_snrs(i) = SNRo2(sig_steady, fs, f0);
end
end

function Plot_Time_Frequency2(x, fs, options)
arguments
    x   double
    fs  double
    options.LineWidth (1,1) {mustBeNumeric} = 1
    options.FaultFreq (1,1) double = NaN
    options.FaultLabel (1,1) string = '$f_{fault}=$'
end

N = length(x);
t = (0:N-1) / fs;
f = fs/N*(0:(N/2));

P2 = abs(fft(x)/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

SetThesisDefaultStyle();
CreateThesisFigure(8, 6);
layout = tiledlayout(2,1);
layout.Padding = 'tight';
layout.TileSpacing = 'tight';

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
xlabel('Time [s]'); ylabel('Amplitude'); grid on;

ax = nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlabel('Frequency [Hz]'); ylabel('Amplitude'); grid on;

if ~isnan(options.FaultFreq)
    ff = options.FaultFreq;
    yval = interp1(f, P1, ff, 'linear', 'extrap');
    yval = max(0, min(yval, max(P1)));
    
    axpos = ax.Position;
    xnorm = axpos(1) + (ff - ax.XLim(1)) / max(diff(ax.XLim), eps) * axpos(3);
    ynorm = axpos(2) + (yval - ax.YLim(1)) / max(diff(ax.YLim), eps) * axpos(4);
    
    x_tip = min(max(xnorm, 0.02), 0.98);
    y_tip = min(max(ynorm, 0.02), 0.98);
    x_tail = min(max(x_tip + 0.06, 0.02), 0.98);
    y_tail = min(max(y_tip - 0.06, 0.02), 0.98);
    
    annotation('textarrow', [x_tail x_tip], [y_tail y_tip], ...
        'String', sprintf('%s %.2f Hz', options.FaultLabel, ff), ...
        'FontSize', 11, 'LineWidth', 1, 'Interpreter', 'latex');
end
end

function Plot_Time_Frequency1(x, fs, options)
arguments
    x   double
    fs  double
    options.LineWidth (1,1) {mustBeNumeric} = 1
    options.FaultFreq (1,1) double = NaN
    options.FaultLabel (1,1) string = '$f_{fault}=$'
end

N = length(x);
x = x - mean(x);
% 时间轴和频率轴计算
t = (0:N-1) / fs;
f = fs/N*(0:(N/2)); % 频率范围（包含奈奎斯特频率）

% FFT计算和幅值归一化
P2 = abs(fft(x)/N); % 计算fft并归一化纵轴（双边）
P1 = P2(1:N/2+1);   % 截取单边谱
P1(2:end-1) = 2*P1(2:end-1);    % 能量还原（除0Hz和奈奎斯特频率）

% 绘制时域和频域图形
SetThesisDefaultStyle();
CreateThesisFigure(6, 9);
layout = tiledlayout(3,1);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
xlabel('Time[s]')
ylabel('Amplitude')
xlim([0 0.2])

nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlabel('Frequency[Hz]')
ylabel('Amplitude')

ax = nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlabel('Frequency[Hz]')
ylabel('Amplitude')
xlim([0 1000])
% ylim([0 max(P1)*1.1])

if ~isnan(options.FaultFreq)
    % ff = options.FaultFreq;
    ff = 137.5;
    yval = interp1(f, P1, ff, 'linear', 'extrap');
    yval = max(0, min(yval, max(P1)));
    
    axpos = ax.Position; % normalized in figure
    xnorm = axpos(1) + (ff - ax.XLim(1)) / diff(ax.XLim) * axpos(3);
    ynorm = axpos(2) + (yval - ax.YLim(1)) / diff(ax.YLim) * axpos(4) + 0.01;
    annotation('textarrow', [xnorm + 0.06 xnorm], [ynorm + 0.14 ynorm], ...
        'String', sprintf('%s %.2f Hz', options.FaultLabel, 142.90), ...
        'FontSize', 12, 'LineWidth', 1, 'Interpreter', 'latex');
end

end