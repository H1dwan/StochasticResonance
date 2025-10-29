function y = CBSR_Potential(x, varargin)
% CBSR_Potential 计算经典双稳态随机共振(CBSR)系统的势函数
% 输入:
%   x - 自变量（可以是标量或向量）
%   a - 参数，默认值 1
%   b - 参数，默认值 1
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'a', 1);
addParameter(p, 'b', 1);
parse(p, varargin{:});

% 提取参数
a = p.Results.a;
b = p.Results.b;

% 计算势函数值 U(x) = -a*x^2/2 + b*x^4/4
y = -a * x.^2 / 2 + b * x.^4 / 4;
end