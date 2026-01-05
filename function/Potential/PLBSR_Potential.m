function U = PLBSR_Potential(x, U0, L0)
% PLBSR_Potential 计算分段线性双稳态随机共振(PLBSR, Piecewise Linear Bistable Stochastic Resonance)系统的势函数
%
% 输入:
%   x - 自变量（可以是标量或向量）
%   U0 - 势阱深度（可选，默认值为1）
%   L0 - 势阱宽度（可选，默认值为1）
%
% 输出: y - 对应的势函数值
%
% 来源：基于分段线性非饱和随机共振的机械早期故障诊断方法研究

if nargin < 2
    U0 = 1; % 默认值
    L0 = 1; % 默认值
end

% 分段计算势函数值
U = zeros(size(x));

% 区间1: x < -L0
idx1 = x < -L0;
U(idx1) = -U0/L0 .* (x(idx1) + L0);

% 区间2: -L0 <= x < 0
idx2 = x >= -L0 & x < 0;
U(idx2) = U0 * (x(idx2)/L0 + 1);

% 区间3: 0 <= x < L0
idx3 = x >= 0 & x < L0;
U(idx3) = U0/L0 * (-x(idx3)/L0 + 1);

% 区间4: x >= L0
idx4 = x >= L0;
U(idx4) = U0/L0 * (x(idx4) - L0);
end