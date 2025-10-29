function dxdt = CBSR_Dynamics(x, a, b)
% CBSR_Dynamics: 经典双稳态系统动力学方程
% dx/dt = a*x - b*x^3 + s_in
%
% 输入参数:
%   x: 当前状态
%   a, b: 系统参数
% 输出参数:
%   dxdt: 状态变化率

dxdt = a*x - b*x^3;

end