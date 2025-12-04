function snr_db  = SNRo(x, fs, f0)
% SNRo  计算单一目标频率分量的输出信噪比（dB）
%
%   snr_db = SNRo(x, fs, f0)
%
% 输入参数：
%   x   : 输出信号向量（行/列均可）
%   fs  : 采样频率 (Hz)
%   f0  : 目标信号中心频率 (Hz)，用于从功率谱中取"信号"功率
%
% 输出参数：
%   snr_db : 输出信噪比，单位 dB
%
% 说明：
%   - 本函数假定输出中主要有一个目标频率分量 f0，其余视为噪声。
%   - 信号功率通过累加目标频率附近几个频率点的功率来估算；
%   - 噪声功率取除这些点及其邻居若干点以外所有谱线的平均功率。

% 输入参数验证
if ~isnumeric(x) || isempty(x)
    error('SNRo:InvalidSignal', '输入信号必须是非空数值向量');
end

if ~isnumeric(fs) || fs <= 0
    error('SNRo:InvalidSamplingRate', '采样频率必须是正数');
end

if ~isnumeric(f0) || f0 <= 0 || f0 > fs/2
    error('SNRo:InvalidTargetFrequency', '目标频率必须在(0, fs/2]范围内');
end

% 去均值，避免 DC 分量过大影响谱估计
x = x - mean(x);

n = length(x);

% 使用汉宁窗减小频谱泄漏（也可改为 boxcar 窗）
% w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));  % Hann window
xw = x .* 1;

% 计算 FFT
x_fft = fft(xw);

% 双边功率谱（幅度平方 / N）
p2 = (abs(x_fft) / n).^2;

% 单边功率谱（0 ~ fs/2）
half_n = floor(n/2);
p1 = p2(1:half_n+1);
% 非 DC 与非 Nyquist 频率处的功率乘 2（能量对称）
if half_n > 1
    p1(2:end-1) = 2 * p1(2:end-1);
end

% 对应的频率轴
f = (0:half_n) * fs / n;

% 找到最接近 f0 的频率索引
[~, idx0] = min(abs(f - f0));

% 信号功率：目标频率点及其相邻点的功率总和
% 使用3个bin（中心bin和左右各一个）来估算信号功率，缓解频谱泄漏影响
signal_bins_start = max(1, idx0-1);
signal_bins_end = min(length(p1), idx0+1);
signal_power = sum(p1(signal_bins_start:signal_bins_end));

% 噪声功率：除信号点及其左右若干个"保护带"之外的平均谱功率
guard_bins = 2;  % 保护带宽度，可视情况修改
noise_mask = true(size(p1));
% 去掉 DC
noise_mask(1) = false;
% 去掉 f0 附近的保护带
noise_mask(max(1, idx0-guard_bins):min(length(p1), idx0+guard_bins)) = false;

% 若都被去掉则给出回退处理
if ~any(noise_mask)
    % 极端情况下，直接将噪声功率设为非常小的值
    noise_power = eps;
else
    noise_power = sum(p1(noise_mask));
end

% 防止除零
if noise_power <= 0
    noise_power = eps;
end

% 计算 SNR（dB）
snr_db = 10 * log10(signal_power / noise_power);

end