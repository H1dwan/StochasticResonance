function h_spec = SpecEntropy(x, ~)
% SpecEntropy  计算信号的归一化谱熵（0~1）
%
%   h_spec = SpecEntropy(x, fs)
%
% 输入参数：
%   x  : 信号向量
%   fs : 采样频率 (Hz)，主要用于频带解释，这里只用来构造频率轴
%
% 输出参数：
%   h_spec : 归一化谱熵，范围约为 [0, 1]
%       - 越接近 1：谱越“平坦”，接近噪声
%       - 越接近 0：谱能量越集中，接近单/少数频率分量
%
% 说明：
%   - 本实现采用简单的 FFT + Hann 窗做功率谱估计。
%   - 若需要更精细的谱估计，可改为使用 pwelch 等。

% 转为列向量
x = x(:);

% 去均值减少 DC 分量
x = x - mean(x);

n = length(x);

% Hann 窗
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
xw = x .* w;

% FFT
x_fft = fft(xw);

% 双边功率谱
p2 = (abs(x_fft) / n).^2;

% 单边功率谱
half_n = floor(n/2);
p1 = p2(1:half_n+1);
if half_n > 1
    p1(2:end-1) = 2 * p1(2:end-1);
end

% 可选：去掉 DC (如果你认为 DC 不属于“有效频带”)
% 这里保留 DC，仅在归一化时一起考虑。

% 对功率谱归一化，得到“频域概率分布”
psd_sum = sum(p1);
if psd_sum <= 0
    h_spec = 0;
    return;
end
p = p1 / psd_sum;

% 为避免 log(0)，将零概率替换为一个极小值
p(p == 0) = eps;

% 未归一化的谱熵（以 log2 为底，单位：bits）
h = -sum(p .* log2(p));

% 最大熵（所有频点概率均等）
h_max = log2(length(p));

% 归一化谱熵
h_spec = h / h_max;
end