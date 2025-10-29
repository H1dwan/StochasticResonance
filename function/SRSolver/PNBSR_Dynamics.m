function dxdt = PNBSR_Dynamics(x, a, b, c)
% MUBSR_Dynamics: 多稳态双稳态系统动力学方程
% 根据给定的分段势函数计算状态变化率
%
% 输入参数:
%   x: 当前状态
%   a, b: 结构参数
%   c: 结构参数，控制开口大小
% 输出参数:
%   dxdt: 状态变化率

% 计算分界点
x0 = sqrt(a/b);

if x < -x0
    dxdt = -c * (x + x0);
elseif x <= x0
    dxdt = a*x - b*x^3;
else    
    dxdt = -c * (x - x0);
end