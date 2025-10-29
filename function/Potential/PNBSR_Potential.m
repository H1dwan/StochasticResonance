function y = PNBSR_Potential(x, varargin)
% PNBSR_potential 计算分段非线性双稳态随机共振(PNBSR, Piecewise Non-Linear Bistable Stochastic Resonance)
% 系统的势函数
%
% 输入:
%   x - 自变量（可以是标量或向量）
%   a3 - 参数，默认值 1
%   b3 - 参数，默认值 1
%   c3 - 参数，默认值 1
%
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'a3', 1);
addParameter(p, 'b3', 1);
addParameter(p, 'c3', 1);
parse(p, varargin{:});

% 提取参数
a3 = p.Results.a3;
b3 = p.Results.b3;
c3 = p.Results.c3;

% 计算关键点
threshold = sqrt(a3 / b3);

% 分段计算势函数值
y = zeros(size(x));
% 区间1: x < -sqrt(a3/b3)
idx1 = x < -threshold;
y(idx1) = 0.5 * c3 * (x(idx1) + sqrt(a3/b3)).^2 - a3^2 / (4*b3);

% 区间2: -sqrt(a3/b3) <= x <= sqrt(a3/b3)
idx2 = abs(x) <= threshold;
y(idx2) = -a3/2 * x(idx2).^2 + b3/4 * x(idx2).^4;

% 区间3: x > sqrt(a3/b3)
idx3 = x > threshold;
y(idx3) = 0.5 * c3 * (x(idx3) - sqrt(a3/b3)).^2 - a3^2 / (4*b3);
end