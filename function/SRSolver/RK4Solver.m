function s_out = RK4Solver(dynamicsFunc, s_in, h, varargin)
% RK4Solver: 通用四阶龙格-库塔法求解器，适用于随机共振系统
% 输入参数:
%   dynamicsFunc: 函数句柄，定义系统的动力学方程，格式：dxdt = dynamicsFunc(x, params...)
%   s_in: 输入信号
%   h: 数值积分的时间步长，通常为采样频率的倒数
%   varargin: 可变参数
%     第一个可选参数可以是初始条件 (默认为0)
%     其余参数会传递给动力学函数
% 输出参数:
%   s_out: 系统输出信号，即x(t)

% 提取初始条件
if ~isempty(varargin)
    x0 = varargin{1};
else
    x0 = 0;
end

N = length(s_in);
s_out = zeros(N, 1);
s_out(1) = x0;

for i = 1 : N-1
    % 计算四阶龙格-库塔法的四个中间值
    k1 = h * dynamicsFunc(s_out(i)) + h * s_in(i);
    k2 = h * dynamicsFunc(s_out(i) + 0.5*k1) + h * s_in(i);
    k3 = h * dynamicsFunc(s_out(i) + 0.5*k2) + h * s_in(i+1);
    k4 = h * dynamicsFunc(s_out(i) + k3) + h * s_in(i+1);
    
    % 更新状态
    s_out(i+1) = s_out(i) + (k1 + 2*k2 + 2*k3 + k4) / 6;

    if abs(s_out(i+1)) > 5
        s_out(i+1) = 5 * sign(s_out(i+1)); % 限制输出范围，防止数值爆炸
    end
end

end
