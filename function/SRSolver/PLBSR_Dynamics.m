function drift = PLBSR_Dynamics(x, U0, L0)
% PLBSR_Dynamics: 分段线性双稳态系统的漂移函数
%
% 输入参数:
%   x: 当前状态
%   a: 结构参数，控制开口大小（x轴截距）
%   b: 结构参数，控制势阱位置xm
%   c: 结构参数，控制势垒高度ΔU
% 输出参数:
%   drift: 漂移项

% 根据x所在区间选择不同的力函数
mask_left_outer = x < -L0;
mask_left_inner = x >= -L0 & x < 0;
mask_right_inner = x >= 0 & x < L0;
mask_right_outer = x >= L0;

drift = zeros(size(x));

drift(mask_left_outer) = U0 / L0;
drift(mask_left_inner) = -U0 / L0;
drift(mask_right_inner) = U0 / L0;
drift(mask_right_outer) = -U0 / L0;

end