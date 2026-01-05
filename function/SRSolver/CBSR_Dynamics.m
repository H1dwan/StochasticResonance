function drift = CBSR_Dynamics(x, a, b)
% CBSR_Dynamics: 经典双稳态系统的漂移函数
% 
% 语法:
%   drift = CBSR_Dynamics(x, a, b)
%
% 输入参数:
%   x: 当前状态
%   a, b: 系统参数
%
% 输出参数:
%   drift: 漂移项

drift = a*x - b*x.^3;
end