function h_sym = SymbEntropy(x, window_length, threshold)
% SymbEntropy  计算二值符号序列的归一化香农熵（0~1）
%
%   h_sym = SymbEntropy(x, window_length, threshold)
%
% 输入参数：
%   x             : 实值信号向量
%   window_length : 窗口长度 L（符号“单词”长度），例如 4 或 5
%                   若省略，则默认 5
%   threshold     : 二值化阈值，默认使用 mean(x)
%
% 输出参数：
%   h_sym : 归一化符号熵，范围约为 [0, 1]
%       - 越接近 1：符号序列越“随机”，模式分布均匀
%       - 越接近 0：序列越“规则”，模式集中
%
% 说明：
%   1. 将 x 二值化为 0/1 序列；
%   2. 用长度为 L 的滑动窗口构造二进制“单词”；
%   3. 每个单词看作一个 L 位二进制数（左高位），映射为整数 0~2^L-1；
%   4. 统计所有状态频率，计算香农熵，并除以最大熵 L 得到归一化结果。

% ===== 参数处理 =====
if nargin < 2 || isempty(window_length)
    window_length = 5;  % 默认窗口长度
end

if nargin < 3 || isempty(threshold)
    threshold = mean(x);  % 默认阈值为均值
end

% 转为列向量
x = x(:);
n = length(x);

if n < window_length
    h_sym = NaN;
    return;
end

% ===== 第一步：二值化 =====
% 大于等于阈值设为 1，小于阈值设为 0
sym_seq = x >= threshold;
sym_seq = uint8(sym_seq);   % 存成 0/1

% ===== 第二步：构造窗口“单词”并映射为整数 =====
num_words = n - window_length + 1;

% 保存每个窗口对应的整数 ID（0 ~ 2^L-1）
word_ids = zeros(num_words, 1, 'uint32');

% 生成二进制权重，如 L=4 -> [8; 4; 2; 1]
% 注意：这里用 double，方便后面做点积
weights = 2.^(window_length-1:-1:0).';   % 列向量 (L×1)

for i = 1:num_words
    % 当前窗口 L×1 的 0/1 向量
    window_bits = sym_seq(i:i+window_length-1);
    
    % 显式转 double 做点积，结果是一个标量
    word_val = double(window_bits).' * weights;  % 1×L · L×1 -> 标量
    
    % 存到 uint32 整数 ID，中间结果理论上在 [0, 2^L-1]
    word_ids(i) = uint32(word_val);
end

% ===== 第三步：统计每个状态出现的频率 =====
num_states = 2^window_length;

counts = zeros(num_states, 1);
% word_ids 范围 [0, num_states-1]，索引时 +1
for i = 1:num_words
    idx = double(word_ids(i)) + 1;
    counts(idx) = counts(idx) + 1;
end

probs = counts / num_words;

% 去掉为 0 的概率，避免 log(0)
probs(probs == 0) = [];

% ===== 第四步：香农熵及归一化 =====
h = -sum(probs .* log2(probs));   % bits
h_max = window_length;            % 最大熵：2^L 等概率 -> L bits

h_sym = h / h_max;
end
