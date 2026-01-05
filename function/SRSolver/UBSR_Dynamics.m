function drift = UBSR_Dynamics(x, a, b)
% UBSR_Dynamics  分段构造的抗饱和双稳势函数的漂移项
%   drift = -dU/dx，其中U为UBSR势函数
%
% 功能说明：
%   计算随机微分方程中的漂移项，用于双稳态系统动力学分析
%
% 输入参数：
%   'x'     - 位置坐标
%   varargin - 可选参数名称-值对，包括:
%       'a' - 正实数参数，控制势阱深度和宽度
%       'b' - 正实数参数，影响势函数的非线性程度
%
% 输出参数：
%   drift - 漂移函数值（与输入x维度相同）
%
% 数学模型：
%   dU/dx = {
%       -a^2/(4b) * 1/(c - sqrt(a/b)), x < -sqrt(a/b)
%       -a*x + b*x^3,                 -sqrt(a/b) <= x <= sqrt(a/b)
%       a^2/(4b) * 1/(c - sqrt(a/b)), x > sqrt(a/b)
%   }

c = sqrt(2*a/b);

% 计算分界点
x_left = -sqrt(a/b);
x_right = sqrt(a/b);

% 根据x所在区间选择不同的力函数（支持向量化）
drift = zeros(size(x));

mask_left = x < x_left;
mask_mid = x >= x_left & x <= x_right;
mask_right = x > x_right;

drift(mask_left) = a^2/(4*b) / (c - sqrt(a/b));
drift(mask_mid) = a.*x(mask_mid) - b.*x(mask_mid).^3;
drift(mask_right) = -a^2/(4*b) / (c - sqrt(a/b));

end