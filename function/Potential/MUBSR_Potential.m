function y = MUBSR_Potential(x, varargin)
% MUBSR_Potential 计算多参数抗饱和双稳随机共振(MUBSR, Multi-parameter Unsaturation Bistable Stochastic Resonance)
% 系统的势函数
%
% 输入:
%   x - 自变量（可以是标量或向量）
%   a1 - 参数，默认值 1
%   b1 - 参数，默认值 1
%   c1 - 参数，默认值 sqrt(a1/b1)+0.5 (需满足c1 > sqrt(a1/b1))
%
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'a1', 1);
addParameter(p, 'b1', 1);
addParameter(p, 'c1', sqrt(1/1) + 0.5); % 默认值确保c1 > sqrt(a1/b1)
parse(p, varargin{:});

% 提取参数
a1 = p.Results.a1;
b1 = p.Results.b1;
c1 = p.Results.c1;

% 验证参数条件：c1 > sqrt(a1/b1)
if c1 <= sqrt(a1/b1)
    error('参数c1必须大于sqrt(a1/b1)');
end

% 计算关键点
threshold = sqrt(a1 / b1);

% 分段计算势函数值
y = zeros(size(x));
% 区间1: x < -sqrt(a1/b1)
idx1 = x < -threshold;
y(idx1) = -a1^2/(4*b1) .* (x(idx1) + c1) ./ (c1 - sqrt(a1/b1));

% 区间2: -sqrt(a1/b1) <= x <= sqrt(a1/b1)
idx2 = abs(x) <= threshold;
y(idx2) = -a1/2 * x(idx2).^2 + b1/4 * x(idx2).^4;

% 区间3: x > sqrt(a1/b1)
idx3 = x > threshold;
y(idx3) = a1^2/(4*b1) .* (x(idx3) - c1) ./ (c1 - sqrt(a1/b1));
end