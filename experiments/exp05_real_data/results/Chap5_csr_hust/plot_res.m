clc; clear; close all;

%% 1) 加载结果文件
script_dir = fileparts(mfilename('fullpath'));
mat_name = 'res_0.5X_I_20Hz_2026_csr_global_best.mat';
mat_path = fullfile(script_dir, mat_name);

if ~exist(mat_path, 'file')
    error('未找到结果文件: %s', mat_path);
end

S = load(mat_path);
if ~isfield(S, 'results')
    error('文件 %s 中未找到变量 results', mat_path);
end
results = S.results;

required_top = {'config', 'best', 'metrics', 'signals'};
for i = 1:numel(required_top)
    if ~isfield(results, required_top{i})
        error('results.%s 缺失，无法可视化。', required_top{i});
    end
end

cfg = results.config;
best = results.best;
met = results.metrics;
sig = results.signals;

fprintf('=== 加载结果成功 ===\n');
if isfield(cfg, 'data_filename')
    fprintf('Data: %s\n', cfg.data_filename);
end
if isfield(cfg, 'f_fault_real')
    fprintf('Fault Frequency: %.4f Hz\n', cfg.f_fault_real);
end
fprintf('Input SNR (SR): %.4f dB\n', met.input_snr_sr);
fprintf('Final Gain vs Input (SR): %.4f dB\n', met.final_gain_vs_input_sr);

%% 2) 优化收敛曲线
if isfield(best, 'best_curve') && ~isempty(best.best_curve)
    figure('Color', 'w', 'Position', [100, 120, 760, 360]);
    plot(1:numel(best.best_curve), -best.best_curve, 'LineWidth', 1.6);
    grid on;
    xlabel('Iteration');
    ylabel('Best SNR (dB)');
    title('Global PSO Convergence (3-layer HSUBSR)');
end

%% 3) 层级 SNR / 增益柱状图
layer_snr = met.layer_snr_sr(:)';
num_layers = numel(layer_snr);

figure('Color', 'w', 'Position', [120, 170, 780, 360]);
subplot(1, 2, 1);
bar(1:num_layers, layer_snr, 0.6);
grid on;
set(gca, 'XTick', 1:num_layers, 'XTickLabel', compose('Layer-%d', 1:num_layers));
ylabel('SNR (dB)');
title('Layer-wise SNR (SR domain)');

subplot(1, 2, 2);
bar(1:num_layers, layer_snr - met.input_snr_sr, 0.6);
grid on;
set(gca, 'XTick', 1:num_layers, 'XTickLabel', compose('Layer-%d', 1:num_layers));
ylabel('Gain vs Input (dB)');
title('Layer-wise Gain (SR domain)');

%% 4) SR 域输入与各层输出时域波形
if isfield(cfg, 'fs_sr') && ~isempty(cfg.fs_sr)
    fs_sr = cfg.fs_sr;
else
    fs_sr = 1;
end

input_sr = sig.input_sr(:);
layer_outputs = sig.layer_outputs_sr;

if size(layer_outputs, 1) ~= numel(input_sr)
    error('signals.layer_outputs_sr 与 signals.input_sr 长度不一致。');
end

t_sr = (0:numel(input_sr)-1)' / fs_sr;

figure('Color', 'w', 'Position', [80, 60, 1000, 760]);
subplot(num_layers + 1, 1, 1);
plot(t_sr, input_sr, 'k', 'LineWidth', 1.0);
grid on; ylabel('Amp');
title(sprintf('Input (SR domain), SNR=%.3f dB', met.input_snr_sr));

for k = 1:num_layers
    subplot(num_layers + 1, 1, k + 1);
    plot(t_sr, layer_outputs(:, k), 'LineWidth', 1.0);
    grid on; ylabel('Amp');
    title(sprintf('Layer-%d Output (SR domain), SNR=%.3f dB', k, layer_snr(k)));
end
xlabel('Time (s)');

%% 5) 最终输出频谱（恢复到原频率）
if isfield(cfg, 'rescale_R') && isfield(cfg, 'fs_sr')
    fs_restore = cfg.fs_sr * cfg.rescale_R;
else
    fs_restore = 1;
end

x_final = sig.final_output_sr(:);
f_fault = NaN;
if isfield(cfg, 'f_fault_real')
    f_fault = cfg.f_fault_real;
end

Plot_Time_Frequency2(x_final, fs_restore, 'FaultFreq', f_fault, 'FaultLabel', '$f=$');

%% 6) 参数打印
if isfield(best, 'best_params')
    fprintf('Best params = %s\n', mat2str(best.best_params, 4));
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
    options.FaultLabel (1,1) string = '$f_{BPFI}=$'
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
CreateThesisFigure(9, 3);
layout = tiledlayout(1,2);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
xlabel('Time[s]')
ylabel('Amplitude')
xlim([0 0.2])

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