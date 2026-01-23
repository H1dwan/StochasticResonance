function [snr_theory_freq, snr_adiabatic, r_k] = HSUBSR_SNR(xm, dU, shape, f, A, D_list)
%HSUBSR_SNR  计算HSUBSR模型的信噪比理论值
%
% 语法:
%   [snr_theory_freq, snr_adiabatic, r_k] = HSUBSR_SNR(xm, dU, shape, f, A, D_list)
%
% 输入参数:
%   xm      - 势阱位置 (标量)
%   dU      - 势垒高度 (标量)
%   shape   - 形状参数，传递给CalibrateHSUBSR (标量)
%   f       - 激励频率 (Hz)
%   A       - 激励幅值 (与模型量纲一致)
%   D_list  - 噪声强度列表 (向量)
%
% 输出参数:
%   snr_theory_freq - 频率修正后的理论信噪比 (与D_list同长度向量)
%   snr_adiabatic   - 绝热近似信噪比 (与D_list同长度向量)
%   r_k             - Kramers跃迁率 (与D_list同长度向量)
%
% 说明:
%   先通过CalibrateHSUBSR得到双稳势参数，计算小扰动振荡角频率w0、势阱角
%   频率wb，再使用Kramers公式求跃迁率r_k。绝热近似信噪比随r_k、噪声强度
%   和驱动幅值计算，频率修正采用Lorentzian Roll-off因子 1 / (1 + (pi*f/r_k)^2)。

[a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape);
wb = sqrt(a*k1 - b*k2);
term_w0 = b*k2*(sech(k2*xm))^2 - a*k1*(sech(k1*xm))^2;
w0 = sqrt(abs(term_w0));

r_k = zeros(length(D_list), 1);
snr_adiabatic = zeros(length(D_list), 1);
snr_theory_freq = zeros(length(D_list), 1); % B. 频率修正理论 (Lorentzian Roll-off)

for i = 1:length(D_list)
    D_fixed = D_list(i);
    r_k(i) = (w0 * wb / (2*pi)) * exp(-dU / D_fixed);
    % snr_adiabatic(i) = (xm^2 * A^2 * r_k(i)) / (4 * D_fixed^2);
    snr_adiabatic(i) = (pi/2) * (A*xm/D_fixed)^2 * r_k(i); % A. 绝热近似理论
    correction = 1 / (1 + (pi * f / r_k(i))^2); % 频率修正因子
    snr_theory_freq(i) = snr_adiabatic(i) * correction;
end

end