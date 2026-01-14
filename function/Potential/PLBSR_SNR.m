function [snr_theory_freq, snr_adiabatic, r_k] = PLBSR_SNR(U0, L0, f, A, D_list)
%PLBSR_SNR 计算分段线性双势阱随机共振(PLBSR)系统的理论信噪比
%
%   PLBSR_SNR 计算周期势双势阱随机共振系统在给定参数下的信噪比，
%   包括绝热近似下的信噪比以及考虑频率修正因子(Lorentzian Roll-off)后的理论值。
%
%   语法:
%       [snr_theory_freq, snr_adiabatic, r_k] = PLBSR_SNR(U0, L0, f, A, D_list)
%
%   输入参数:
%       U0          - 势垒高度，标量，表示两个势阱之间的能量差
%       L0          - 势阱位置 xm
%       f           - 输入信号频率 (Hz)，标量
%       A           - 输入信号振幅，标量
%       D_list      - 噪声强度列表 (列向量)，包含待计算的扩散系数值
%
%   输出参数:
%       snr_theory_freq - 应用频率修正后的理论信噪比 (行向量)
%       snr_adiabatic   - 绝热近似下的信噪比 (行向量)
%       r_k             - Kramers逃逸速率 (行向量)，按噪声强度D_list对应
%
%   说明:
%       本函数基于随机共振理论，计算周期双势阱系统在弱周期驱动和
%       白噪声作用下的信噪比。其中：
%       - r_k：Kramers逃逸速率，表征粒子在两势阱间跳跃的平均速率
%       - snr_adiabatic：绝热极限下的信噪比
%       - snr_theory_freq：考虑系统有限频率响应的修正信噪比
%
%   参考文献:
%       关于周期双势阱随机共振的理论基础，见相关文献
%
%   另见: UBSR_SNR, HSUBSR_SNR

prefactor = U0^2 / (2*L0^2);

r_k = zeros(length(D_list), 1);
snr_adiabatic = zeros(length(D_list), 1);
snr_theory_freq = zeros(length(D_list), 1); % B. 频率修正理论 (Lorentzian Roll-off)

for i = 1:length(D_list)
    D_fixed = D_list(i);
    r_k(i) = prefactor / D_fixed * exp(-U0 / D_fixed);
    
    snr_adiabatic(i) = (pi*A^2*U0^2) / (2*D_fixed^3) * exp(-U0 / D_fixed);
    correction = 1 / (1 + (pi * f / r_k(i))^2); % 频率修正因子
    snr_theory_freq(i) = snr_adiabatic(i) * correction;
end

end