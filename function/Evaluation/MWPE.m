function [mwpe_index, wpe_vec, scales] = MWPE(y, m, tau, max_scale)
% calculate_MWPE 计算多尺度加权排列熵 (Multi-scale Weighted Permutation Entropy)
% 输入:
%   y         - 输入信号
%   m         - 嵌入维数 (推荐 3~5)
%   tau       - 延迟时间 (通常设为 1，因为粗粒化已经处理了时间尺度)
%   max_scale - 最大尺度因子 (推荐 20~50，取决于 fs/f0)
% 输出:
%   mwpe_index - 综合指标 (通常取所有尺度的平均值，或特定频带的平均值)
%   wpe_vec    - 各个尺度下的 WPE 值向量
%   scales     - 尺度向量 1:max_scale

if nargin < 4, max_scale = 20; end
if nargin < 3, tau = 1; end
if nargin < 2, m = 3; end

wpe_vec = zeros(max_scale, 1);
scales = 1:max_scale;

for s = 1:max_scale
    % 1. 粗粒化处理
    y_s = CoarseGraining(y, s);
    
    % 2. 计算当前尺度下的 WPE
    % 注意：随着 s 增大，序列变短，需确保长度足够
    if length(y_s) < 5 * factorial(m)
        warning('尺度 %d 下序列过短，可能导致熵值估计不准', s);
        wpe_vec(s) = NaN;
    else
        wpe_vec(s) = local_WPE(y_s, m, tau);
    end
end

% 3. 计算综合指标 (Innovation Point)
% 策略：取中大尺度 (Scale 5 到 max_scale) 的平均值
% 理由：小尺度 (1-4) 仍包含较多高频噪声，大尺度主要反映宏观周期结构
valid_indices = 5:max_scale;
if max_scale < 5, valid_indices = 1:max_scale; end

% 忽略 NaN 值
valid_wpe = wpe_vec(valid_indices);
mwpe_index = mean(valid_wpe(~isnan(valid_wpe)));
end

% --- 内部辅助函数: WPE 计算 ---
function WPE = local_WPE(y, m, tau)
% 标准 WPE 计算 (精简版)
y = y(:);
N = length(y);
n_vec = N - (m-1)*tau;

% 重构
X = zeros(n_vec, m);
for i = 1:m
    X(:, i) = y(1+(i-1)*tau : n_vec+(i-1)*tau);
end

% 权重 (方差)
weights = var(X, 0, 2);

% 模式识别
[~, ord] = sort(X, 2);
[~, ~, pattern_ids] = unique(ord, 'rows');

% 加权概率
num_patterns = factorial(m);
p_w = zeros(num_patterns, 1);
total_weight = sum(weights);

if total_weight == 0, WPE = 0; return; end

unique_ids = unique(pattern_ids);
for k = 1:length(unique_ids)
    id = unique_ids(k);
    p_w(k) = sum(weights(pattern_ids == id));
end
p_w = p_w / total_weight;

% 熵
p_w = p_w(p_w > 0);
H_w = -sum(p_w .* log(p_w));
WPE = H_w / log(factorial(m));
end