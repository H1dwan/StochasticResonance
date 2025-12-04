function U = HESUBSR_Potential(x, varargin)
% HESUBSR_Potential  超曲率指数型平滑不饱和双稳势函数
%   U(x) = α[ ln cosh(β(x-d)) + ln cosh(β(x+d)) ] - γ ln cosh(βx)
%
% 详细说明：
%   该函数实现超曲率指数型平滑不饱和双稳势模型，用于随机共振研究。
%
% 输入参数:
%   x - 位置坐标（标量或向量）
%   varargin - 可选参数名称-值对，包括:
%       'alpha' - 势垒高度尺度，α 越大两侧的势阱越深
%       'beta'  - 超曲率参数，增大 β 势阱更窄、局部曲率越大；减小 β 势阱更宽更浅
%       'd'     - 稳态点间距控制量，调节两势阱之间的距离
%       'gamma' - 中心势垒 + 势壁斜率的耦合参数，γ 越大 → 中心势垒高，外侧势壁斜率
%
% 输出参数:
%   U - 计算得到的势能值，与输入x维度相同
%
% 示例:
%   U = HESUBSR_Potential(0, 'alpha', 2, 'd', 0.5);

p = inputParser;
addParameter(p, 'alpha', 1);
addParameter(p, 'beta',  1);
addParameter(p, 'd',     1);
addParameter(p, 'gamma', 1.5);
parse(p, varargin{:});

alpha = p.Results.alpha;
beta  = p.Results.beta;
d     = p.Results.d;
gamma = p.Results.gamma;

U = alpha * ( log(cosh(beta * (x - d))) + log(cosh(beta * (x + d))) ) ...
    - gamma * log(cosh(beta * x));
end