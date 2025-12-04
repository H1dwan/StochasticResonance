function drift = HESUBSR_Dynamics(x, alpha, beta, d, gamma)
% HESUBSR_Dynamics  超曲率指数型平滑不饱和双稳势模型的漂移函数
%   drift = -dU/dx，其中U为HESUBSR势函数
%
% 功能说明：
%   计算随机微分方程中的漂移项，用于双稳态系统动力学分析
%
% 输入参数：
%   'x'     - 位置坐标（标量或向量）
%   'alpha' - 势垒高度尺度，α 越大两侧的势阱越深
%   'beta'  - 超曲率参数，增大 β 势阱更窄、局部曲率越大；减小 β 势阱更宽更浅
%   'd'     - 稳态点间距控制量，调节两势阱之间的距离
%   'gamma' - 中心势垒 + 势壁斜率的耦合参数，γ 越大 → 中心势垒高，外侧势壁斜率
%
% 输出参数：
%   drift - 漂移函数值（与输入x维度相同）
%
% 数学模型：
%   U' = αβ[tanh(β(x-d))+tanh(β(x+d))] - γβ tanh(βx)
%
% 示例：
%   x = -5:0.1:5;
%   drift = HESUBSR_Drift(x, 2, 1.5, 0.8, 1.2);

drift = -alpha * beta * (tanh(beta*(x - d)) + tanh(beta*(x + d))) ...
        + gamma * beta * tanh(beta * x);
end