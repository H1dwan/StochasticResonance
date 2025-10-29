function y = PLBSR_Potential(x, varargin)
% PLBSR_Potential 计算分段线性双稳态随机共振(PLBSR, Piecewise Linear Bistable Stochastic Resonance)系统的势函数
%
% 输入:
%   x - 自变量（可以是标量或向量）
%   c2 - 参数，默认值 1
%   b2 - 参数，默认值 1
%   a2 - 参数，默认值 1
%
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'c2', 1);
addParameter(p, 'b2', 1);
addParameter(p, 'a2', 1);
parse(p, varargin{:});

% 提取参数
c2 = p.Results.c2;
b2 = p.Results.b2;
a2 = p.Results.a2;

% 计算关键点
threshold = b2;

% 分段计算势函数值
y = zeros(size(x));
% 区间1: x < -b2
idx1 = x < -threshold;
y(idx1) = -c2/(a2 - b2) .* (x(idx1) + a2);

% 区间2: -b2 <= x < 0
idx2 = x >= -threshold & x < 0;
y(idx2) = c2/b2 * x(idx2);

% 区间3: 0 <= x < b2
idx3 = x >= 0 & x < threshold;
y(idx3) = -c2/b2 * x(idx3);

% 区间4: x >= b2
idx4 = x >= threshold;
y(idx4) = c2/(a2 - b2) .* (x(idx4) - a2);
end