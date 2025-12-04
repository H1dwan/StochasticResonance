function [alpha, gamma, xw, DeltaU] = Calibrate_HESUBSR(beta, d, rho, DeltaU_target)
% HESUBSR模型参数标定函数
% 该函数通过给定的目标势垒高度数值来标定HESUBSR模型的alpha和gamma参数
%
% 输入参数:
%   beta           - 超曲率参数 β，可取1
%   d              - 势函数参数d，决定双稳态位置
%   rho            - ρ = γ/α，需满足 2*sech^2(βd) < ρ < 2，通常可取1.5
%   DeltaU_target  - 目标势垒高度
%
% 输出参数:
%   alpha          - 标定得到的α参数
%   gamma          - 标定得到的γ参数，满足关系 γ = ρ * α
%   xw             - 右侧势阱位置
%   DeltaU         - 实际势垒高度，等于目标势垒高度DeltaU_target
%
% 势函数形式（单位alpha=1）：
%   U1(x) = ln cosh(β(x-d)) + ln cosh(β(x+d)) - ρ ln cosh(βx)
% 给定 α, γ = ρ α 后：
%   U(x) = α * U1(x)

% ---- 1. 求右侧势阱的位置 xw (>0) ----
% 驻点条件 U'(x)=0 => F(x)=0
F = @(x) tanh(beta*(x-d)) + tanh(beta*(x+d)) - rho * tanh(beta*x);

% 以 x0 = d 为初值使用 fzero
x0 = d;  % 初始猜测
options = optimset('Display','off');
xw = fzero(F, x0, options);   % 右侧稳态点

% ---- 2. 在 α=1 情况下计算势垒高度 ΔU_unit ----
U1_0 = 2 * log(cosh(beta*d)); % x=0
U1_w = log(cosh(beta*(xw-d))) + log(cosh(beta*(xw+d))) ...
    - rho * log(cosh(beta*xw));

DeltaU_unit = U1_0 - U1_w;    % 对应α=1时的势垒高度

% ---- 3. 令 α 使得 ΔU = ΔU_target ----
alpha = DeltaU_target / DeltaU_unit;
gamma = rho * alpha;
DeltaU = DeltaU_target;  % 标定后实际势垒高度

end