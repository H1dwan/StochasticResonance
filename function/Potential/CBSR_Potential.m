function y = CBSR_Potential(x, a, b)
% CBSR_Potential 计算经典双稳态随机共振(CBSR)系统的势函数
% 输入:
%   x - 自变量（可以是标量或向量）
%   a - 参数，默认值 1
%   b - 参数，默认值 1
% 输出: y - 对应的势函数值

if nargin < 2
    a = 1;
end

if nargin < 3
    b = 1;
end

% 计算势函数值 U(x) = -a*x^2/2 + b*x^4/4
y = -a * x.^2 / 2 + b * x.^4 / 4;
end