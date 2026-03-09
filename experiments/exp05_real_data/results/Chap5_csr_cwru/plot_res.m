clc; clear; close all;

%% 1. 加载结果文件 ==========================================================
res_path = 'res_135_2026_csr_global_58.mat';

S = load(res_path);
if ~isfield(S, 'results')
    error('结果文件中不存在变量 results: %s', res_path);
end
results = S.results;

fprintf('Loaded result: %s\n', res_path);

required_fields = {'signals', 'config', 'metrics'};
for i = 1:numel(required_fields)
    if ~isfield(results, required_fields{i})
        error('results 缺少字段: %s', required_fields{i});
    end
end

if ~isfield(results.signals, 'input_sr') || ~isfield(results.signals, 'layer_outputs_sr')
    error('results.signals 缺少 input_sr 或 layer_outputs_sr 字段，请先重新运行 chap5_csr_cwru.m');
end

%% 2. 数据提取与指标打印 ====================================================
x_in = results.signals.input_sr(:);
X_layers = results.signals.layer_outputs_sr;

if size(X_layers, 1) ~= numel(x_in)
    if size(X_layers, 2) == numel(x_in)
        X_layers = X_layers.';
    else
        error('layer_outputs_sr 尺寸与 input_sr 不一致');
    end
end

n_layers = size(X_layers, 2);
fs_sr = results.config.fs_sr;
R = results.config.rescale_R;
fs_restore = fs_sr * R;

if isfield(results.config, 'f_target_sr')
    f_target_sr = results.config.f_target_sr;
else
    f_target_sr = NaN;
end

if isfield(results.config, 'f_fault_real')
    f_fault_real = results.config.f_fault_real;
else
    f_fault_real = f_target_sr * R;
end

if isfield(results.metrics, 'input_snr_sr')
    snr_in_sr = results.metrics.input_snr_sr;
else
    snr_in_sr = CalcSteadySNR(x_in, fs_sr, f_target_sr);
end

if isfield(results.metrics, 'layer_snr_sr')
    snr_layers_sr = results.metrics.layer_snr_sr;
else
    snr_layers_sr = zeros(1, n_layers);
    for k = 1:n_layers
        snr_layers_sr(k) = CalcSteadySNR(X_layers(:, k), fs_sr, f_target_sr);
    end
end

fprintf('---------------- CSR 输出可视化 ----------------\n');
fprintf('数据文件: %s\n', results.config.data_filename);
fprintf('输入 SNR (SR域): %.4f dB\n', snr_in_sr);
for k = 1:n_layers
    fprintf('第%d层 SNR (SR域): %.4f dB\n', k, snr_layers_sr(k));
end
fprintf('最终增益: %.4f dB\n', snr_layers_sr(end) - snr_in_sr);
fprintf('------------------------------------------------\n');

%% 3. 时域可视化（输入 + 各层输出） ========================================
N = numel(x_in);
t_sr = (0:N-1)' / fs_sr;

show_ratio = 1.0;
show_len = max(10, round(show_ratio * N));
idx_show = 1:show_len;

figure('Color', 'w', 'Position', [100, 60, 1100, 820]);
subplot(n_layers+1, 1, 1);
plot(t_sr(idx_show), x_in(idx_show), 'k', 'LineWidth', 1.0);
grid on;
ylabel('Amp');
title(sprintf('Input (SR domain), SNR=%.3f dB', snr_in_sr));

for k = 1:n_layers
    subplot(n_layers+1, 1, k+1);
    plot(t_sr(idx_show), X_layers(idx_show, k), 'LineWidth', 1.0);
    grid on;
    ylabel('Amp');
    title(sprintf('Layer-%d Output (SR domain), SNR=%.3f dB', k, snr_layers_sr(k)));
end
xlabel('Time (s)');

%% 4. 频域可视化（输入 + 各层输出） ========================================
figure('Color', 'w', 'Position', [120, 80, 1100, 840]);

[f_sr, p_in] = SingleSideSpectrum(x_in, fs_sr);
subplot(n_layers+1, 1, 1);
plot(f_sr, p_in, 'k', 'LineWidth', 1.0); hold on;
if ~isnan(f_target_sr)
    xline(f_target_sr, 'r--', 'LineWidth', 1.2);
end
grid on;
xlim([0, min(0.2, fs_sr/2)]);
ylabel('|X(f)|');
title('Input Spectrum (SR domain)');

for k = 1:n_layers
    [f_sr, p_k] = SingleSideSpectrum(X_layers(:, k), fs_sr);
    subplot(n_layers+1, 1, k+1);
    plot(f_sr, p_k, 'LineWidth', 1.0); hold on;
    if ~isnan(f_target_sr)
        xline(f_target_sr, 'r--', 'LineWidth', 1.2);
    end
    grid on;
    xlim([0, min(0.2, fs_sr/2)]);
    ylabel('|X(f)|');
    title(sprintf('Layer-%d Spectrum (SR domain)', k));
end
xlabel('Frequency (Hz)');

%% 5. 层级增益可视化 =======================================================
gain_vs_input = snr_layers_sr - snr_in_sr;

figure('Color', 'w', 'Position', [160, 120, 760, 360]);
bar(1:n_layers, gain_vs_input, 0.6, 'FaceColor', [0.3, 0.55, 0.85]);
grid on;
set(gca, 'XTick', 1:n_layers, 'XTickLabel', compose('Layer-%d', 1:n_layers));
ylabel('SNR Gain vs Input (dB)');
title('CWRU Cascaded SR: Layer-wise Gain (SR domain)');

%% 6. 最终输出映射回原频率域展示 ===========================================
x_final = X_layers(:, end);
x_final = x_final - mean(x_final);
x_final = 1 * x_final(0.1*end:end);
f_fault_restore = f_target_sr * R;
Plot_Time_Frequency2(x_final, fs_restore, 'FaultFreq', f_fault_restore, 'FaultLabel', '$f_{BPFO}=$');


%% ======================== 辅助函数 =======================================
function snr_val = CalcSteadySNR(x, fs, f0)
x = x(:);
idx0 = round(0.1 * numel(x));
x_st = x(idx0+1:end);
if isnan(f0)
    snr_val = NaN;
else
    snr_val = SNRo2(x_st, fs, f0);
end
end

function [f, p1] = SingleSideSpectrum(x, fs)
x = x(:);
N = numel(x);
X = fft(x);
p2 = abs(X / N);
p1 = p2(1:floor(N/2)+1);
if numel(p1) > 2
    p1(2:end-1) = 2 * p1(2:end-1);
end
f = fs * (0:floor(N/2)) / N;
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