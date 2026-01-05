function potential = HSUBSR_Potential(x, a, b, k1, k2)
% HSUBSR_Potential 双曲光滑抗饱和双稳势函数
%
% 语法：
%   potential = HSUBSR_Potential(x, a, b, k1, k2)
%
% 输入参数：
%   x        - 位置坐标 标量或向量 (double)
%              系统状态变量，表示粒子或振子的位置
%   a        - 第一项的幅度系数 (double, > 0)
%              强度参数，控制排斥项的幅值
%   b        - 第二项的幅度系数 (double, > 0)
%              强度参数，控制吸引项的幅值
%   k1       - 第一项的斜率/陡峭度参数 (double, > 0)
%              尺度参数，控制排斥势的作用范围与陡
%              峭程度
%   k2       - 第二项的斜率/陡峭度参数 (double, > 0)
%              尺度参数，控制吸引势的作用范围与陡峭程度
%
% 输出参数：
%   potential - 计算得到的势能值 (double)

potential = (b / k2) * log(cosh(k2 * x)) - (a / k1) * log(cosh(k1 * x));
end 