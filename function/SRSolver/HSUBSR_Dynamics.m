function drift = HSUBSR_Dynamics(x, a, b, k1, k2)
% HSUBSR_Dynamics 双曲光滑抗饱和双稳势模型的漂移函数
%
% 语法：
%   drift = HSUBSR_Dynamics(x, a, b, k1, k2)
%
% 输入参数：
%   x      - 当前系统状态/位置 (double)
%            系统状态变量，表示粒子或振子的位置
%
%   a      - 第一项的幅度系数 (double, > 0)
%            强度参数，控制排斥项的幅值
%
%   b      - 第二项的幅度系数 (double, > 0)
%            强度参数，控制吸引项的幅值
%
%   k1     - 第一项的斜率/陡峭度参数 (double, > 0)
%            尺度参数，控制排斥势的作用范围与陡峭程度
%
%   k2     - 第二项的斜率/陡峭度参数 (double, > 0)
%            尺度参数，控制吸引势的作用范围与陡峭程度
%
% 输出参数：
%   drift  - 确定性漂移项/势函数导数的近似值 (double)
%            系统的确定性漂移速率，为势函数关于x的负梯度
%
% 注记：
%  双稳态条件：a < b 且 a*k1 > b*k2

drift = a * tanh(k1 * x) - b * tanh(k2 * x);
end