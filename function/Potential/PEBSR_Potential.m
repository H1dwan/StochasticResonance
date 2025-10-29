function y = PEBSR_Potential(x, varargin)
% PEBSR_potential 计算分段指数双稳态随机共振(PEBSR, Piecewise Exponential Bistable Stochastic Resonance)
% 系统的势函数
%
% 输入:
%   x - 自变量（可以是标量或向量）
%   l4 - 参数，默认值 1
%   c4 - 参数，默认值 1
%   a4 - 参数，默认值 1
%   b4 - 参数，默认值 1
%
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'l4', 1);
addParameter(p, 'c4', 1);
addParameter(p, 'a4', 1);
addParameter(p, 'b4', 1);
parse(p, varargin{:});

% 提取参数
l4 = p.Results.l4;
c4 = p.Results.c4;
a4 = p.Results.a4;
b4 = p.Results.b4;

% 计算关键点：sqrt(a4/b4)
threshold = sqrt(a4 / b4);

% 计算 k4 常数项
k4 = 2 * (l4 + a4^2 / (4*b4));

% 分段计算势函数值
y = zeros(size(x));
% 区间1: x < -sqrt(a4/b4)
idx1 = x < -threshold;
y(idx1) = l4 * exp(-c4 * (x(idx1) + sqrt(a4/b4)) / 2) - k4 / 2;

% 区间2: -sqrt(a4/b4) <= x <= sqrt(a4/b4)
idx2 = abs(x) <= threshold;
y(idx2) = -a4/2 * x(idx2).^2 + b4/4 * x(idx2).^4;

% 区间3: x > sqrt(a4/b4)
idx3 = x > threshold;
y(idx3) = l4 * exp(c4 * (x(idx3) - sqrt(a4/b4)) / 2) - k4 / 2;
end