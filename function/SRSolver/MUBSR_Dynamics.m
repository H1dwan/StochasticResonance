function dxdt = MUBSR_Dynamics(x, a, b, c)
% MUBSR_Dynamics: 多稳态双稳态系统动力学方程
% 根据给定的分段势函数计算状态变化率
%
% 输入参数:
%   x: 当前状态
%   a, b: 结构参数
%   c: 结构参数，控制开口大小
% 输出参数:
%   dxdt: 状态变化率

% 验证参数条件：c > sqrt(a/b)
if c <= sqrt(a/b)
    error('参数c必须大于sqrt(a/b)');
end

% 计算分界点
x_left = -sqrt(a/b);
x_right = sqrt(a/b);

% 根据x所在区间选择不同的力函数
if x < x_left
    % 左侧区域: dx/dt = -dU/dx = a/(4*b) * (x + c1)/(c1 - sqrt(a/b))    
    dxdt = a^2/(4*b) / (c - sqrt(a/b));
elseif x <= x_right
    % 中间区域: dx/dt = -dU/dx = a*x - b*x^3
    dxdt = a*x - b*x^3;
else
    % 右侧区域: dxdt = -dU/dx = -a/(4*b) * (x - c1)/(c1 - sqrt(a/b))
    dxdt = -a^2/(4*b) / (c - sqrt(a/b));
end
end