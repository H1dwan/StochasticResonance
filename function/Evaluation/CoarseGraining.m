function y_coarse = CoarseGraining(y, s)
% CoarseGraining 对时间序列进行粗粒化处理
% 输入:
%   y - 输入原始信号 (向量)
%   s - 尺度因子 (整数, s >= 1)
% 输出:
%   y_coarse - 粗粒化后的信号
%
% 原理: 将序列分成长度为 s 的非重叠窗口，计算每个窗口的均值。
% y_j^s = (1/s) * sum(y(i))  for i from (j-1)s+1 to j*s

% 参数检查
if s < 1
    error('尺度因子 s 必须是大于等于 1 的整数');
end

if s == 1
    y_coarse = y;
    return;
end

y = y(:); % 确保列向量
N = length(y);

% 计算粗粒化后的长度
L = floor(N / s);

% 截断多余的数据点
y_cut = y(1 : L*s);

% 重塑矩阵：每列代表一个窗口，行数为 s
% reshape 会按列填充，所以我们需要先转置或谨慎操作
% 这里构建矩阵 M，大小为 [s, L]
M = reshape(y_cut, s, L);

% 计算每列的均值 (即每个窗口的均值)
y_coarse = mean(M, 1)'; % 结果转为列向量
end