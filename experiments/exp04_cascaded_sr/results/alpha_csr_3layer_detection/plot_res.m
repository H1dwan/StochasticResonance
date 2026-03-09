clc; clear; close all;

%% 1) 加载结果文件 =========================================================
result_files = dir('res_alpha_csr_3layer_detection_1.mat');
if isempty(result_files)
    error('当前目录未找到结果文件: res_alpha_csr_3layer_detection*.mat');
end

[~, idx_latest] = max([result_files.datenum]);
data_file = result_files(idx_latest).name;
load(data_file, 'results');

fprintf('加载结果文件: %s\n', data_file);

if ~isfield(results, 'signals') || ~isfield(results.signals, 'layer_output_all')
    error(['结果文件中缺少每层输出信号。请先重新运行脚本: ', ...
        'experiments/exp04_cascaded_sr/scripts/alpha_noise_csr_3layer_detection.m']);
end

%% 2) 选择试验样本 =========================================================
num_trials = results.config.num_trials;
num_layers = results.config.num_layers;
fs = results.config.fs;
f0 = results.config.f0;

% 默认选择三层输出SNR最高的一次试验
[~, trial_id] = max(results.metrics.snr_3layer);
trial_id = max(1, min(num_trials, trial_id));

clean_sig = results.signals.clean_sig;
noisy_input = results.signals.input_noisy_all(:, trial_id);
layer_outputs = results.signals.layer_output_all(:, :, trial_id);

snr_in = SNRo2(noisy_input, fs, f0);
Plot_Time_Frequency1(clean_sig+noisy_input, fs);

snr_layer = zeros(1, num_layers);
for k = 1:num_layers
    xk = layer_outputs(:, k);
    idx0 = round(0.1 * length(xk));
    snr_layer(k) = SNRo2(xk(idx0+1:end), fs, f0);
end

fprintf('选中试验: #%d/%d\n', trial_id, num_trials);
fprintf('SNR(input) = %.4f dB\n', snr_in);
for k = 1:num_layers
    fprintf('SNR(layer-%d) = %.4f dB\n', k, snr_layer(k));
end

%% 3) 每层输出时频图 =======================================================
N = length(noisy_input);
t = (0:N-1)' / fs;
show_len = min(N, round(400 / f0));
idx_show = 1:show_len;

for k = 1:num_layers
    xk = layer_outputs(:, k);
    Plot_Time_Frequency2(xk(idx_show), fs);
    % sgtitle(sprintf('Trial #%d: Layer-%d Output', trial_id, k));
end


%% 5) 增益柱状图（按层） ===================================================
gain_vs_input = snr_layer - snr_in;
figure('Color', 'w', 'Position', [150, 150, 700, 360]);
bar(1:num_layers, gain_vs_input, 0.6);
grid on;
set(gca, 'XTick', 1:num_layers, 'XTickLabel', compose('Layer-%d', 1:num_layers));
ylabel('SNR Gain vs Input (dB)');
title(sprintf('Trial #%d: Layer-wise Detection Gain', trial_id));


%% 辅助函数 ================================================================
function P1 = SingleSideAmp(x, nfft)
x = x(:);
X = fft(x, nfft);
P2 = abs(X / nfft);
P1 = P2(1:floor(nfft/2)+1);
if numel(P1) > 2
    P1(2:end-1) = 2 * P1(2:end-1);
end
end

function Plot_Time_Frequency1(x, fs, options)
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
% ylim([-20 20])
xlabel('Time[s]')
ylabel('Amplitude')
% title('Time Domain')
% set(gca,'FontSize',14,'FontName','Times New Roman');
% xlim([0 N])
xticks(0:200:1000);

nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlim([0 2.5])
xlabel('Frequency[Hz]')
ylabel('Amplitude')
% title('Frequency Domain')
% set(gca,'FontSize',14,'FontName','Times New Roman');

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
xticks(0:200:1000);

ax = nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlim([0 0.3])
xlabel('Frequency[Hz]')
ylabel('Amplitude')

ff = 0.01;
yval = interp1(f, P1, ff, 'linear', 'extrap');
yval = max(0, min(yval, max(P1)));

axpos = ax.Position; % normalized in figure
xnorm = axpos(1) + (ff - ax.XLim(1)) / diff(ax.XLim) * axpos(3);
ynorm = axpos(2) + (yval - ax.YLim(1)) / diff(ax.YLim) * axpos(4) + 0.01;
annotation('textarrow', [xnorm + 0.06 xnorm], [ynorm - 0.06 ynorm], ...
    'String', sprintf('$f=%.2f\\,\\mathrm{Hz}$', ff), ...
    'FontSize', 12, 'LineWidth', 1, 'Interpreter', 'latex');

end

