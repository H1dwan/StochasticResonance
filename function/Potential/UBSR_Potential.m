function U = UBSR_Potential(x, a, b)
% UBSR_Potential  分段构造的抗饱和双稳势函数（UBSR）
% 来源：An adaptive unsaturated bistable stochastic resonance method and its application
% in mechanical fault diagnosis
%   U(x) = {
%       -a^2/(4b) * (x + c)/(c - sqrt(a/b)), x < -sqrt(a/b)
%       -a/2 * x^2 + b/4 * x^4,             -sqrt(a/b) <= x <= sqrt(a/b)
%       a^2/(4b) * (x - c)/(c - sqrt(a/b)), x > sqrt(a/b)
%   }
%
% 详细说明：
%   该函数实现一种基于分段构造的双稳势模型，用于随机共振研究。
%   势函数在区间交界处连续但导数不连续，可能导致系统响应出现直流偏移或失真。
%
% 输入参数:
%   x - 位置坐标（标量或向量）
%   a - 正实数参数，控制势阱深度和宽度
%   b - 正实数参数，影响势函数的非线性程度
%
% 输出参数:
%   U - 计算得到的势能值，与输入x维度相同
%
% 示例:
%   U = UBSR_Potential(x, 1, 1);

if nargin < 2
    a = 1;
end

if nargin < 3
    b = 1;
end

c = sqrt(2*a/b);

% 计算边界点
sqrt_ab = sqrt(a/b);

% 定义三个区间的条件
cond1 = x < -sqrt_ab;           % 区间1: x < -sqrt(a/b)
cond2 = (-sqrt_ab <= x) & (x <= sqrt_ab); % 区间2: -sqrt(a/b) <= x <= sqrt(a/b)
cond3 = x > sqrt_ab;            % 区间3: x > sqrt(a/b)

% 初始化U为全零数组
U = zeros(size(x));

% 区间1: x < -sqrt(a/b)
U(cond1) = -(a^2)/(4*b) .* (x(cond1) + c) ./ (c - sqrt_ab);

% 区间2: -sqrt(a/b) <= x <= sqrt(a/b)
U(cond2) = -(a/2) .* x(cond2).^2 + (b/4) .* x(cond2).^4;

% 区间3: x > sqrt(a/b)
U(cond3) = (a^2)/(4*b) .* (x(cond3) - c) ./ (c - sqrt_ab);
end