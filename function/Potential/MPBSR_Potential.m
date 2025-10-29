function y = MPBSR_Potential(x, varargin)
% MPBSR_Potential 计算分段构造的非线性势函数（MPBSR）
% 输入:
%   x - 自变量（可以是标量或向量）
%   p - 参数，默认值 1
%   c - 参数，默认值 1
%   u - 参数，默认值 1
% 输出: y - 对应的势函数值

% 使用 inputParser 解析输入参数
p = inputParser;
addParameter(p, 'p', 1);
addParameter(p, 'c', 1);
addParameter(p, 'u', 1);
parse(p, varargin{:});

% 提取参数
p_val = p.Results.p;
c_val = p.Results.c;
u_val = p.Results.u;

% 计算边界值
sqrt_u = sqrt(u_val);

% 初始化输出数组
y = zeros(size(x));

% 分段计算势函数值
% 区间 1: x < -sqrt(u)
idx1 = x < -sqrt_u;
y(idx1) = -(p_val / (4 * (c_val - sqrt_u))) * (x(idx1) + c_val);

% 区间 2: -sqrt(u) <= x <= sqrt(u)
idx2 = (x >= -sqrt_u) & (x <= sqrt_u);
y(idx2) = -(p_val / (2 * u_val)) * x(idx2).^2 + (p_val / (4 * u_val^2)) * x(idx2).^4;

% 区间 3: x > sqrt(u)
idx3 = x > sqrt_u;
y(idx3) = (p_val / (4 * (c_val - sqrt_u))) * (x(idx3) - c_val);
end