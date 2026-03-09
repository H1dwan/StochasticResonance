clc; clear; close all;

D_list = 0:0.01:0.5;

a1 = 1; b1 = 1;
[~, snr_adiabatic1, ~] = CBSR_SNR_Theory(a1, b1, 0.1, 0.1, D_list);

a2 = 1; b2 = 0.5;
[~, snr_adiabatic2, ~] = CBSR_SNR_Theory(a2, b2, 0.1, 0.1, D_list);

a3 = 1; b3 = 2;
[~, snr_adiabatic3, ~] = CBSR_SNR_Theory(a3, b3, 0.1, 0.1, D_list);


SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');

CreateThesisFigure();
tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');

ax = nexttile; hold on;
ax.Position = [0.01 0.01 0.98 0.98];

plot(D_list, snr_adiabatic1, 'o-', 'DisplayName', '$a=1, b=1$');
plot(D_list, snr_adiabatic2, 's-', 'DisplayName', '$a=1, b=0.5$');
plot(D_list, snr_adiabatic3, 'd-', 'DisplayName', '$a=0.5, b=2$');
xticks(0:0.1:0.5);

legend('Location', 'northeast');
xlabel('$D$');
ylabel('SNR');



function [snr_theory_freq, snr_adiabatic, r_k] = CBSR_SNR_Theory(a, b, f, A, D_list)
% CBSR_SNR_Theory 计算经典双稳态随机共振(CBSR)的理论信噪比
%
% 参考文献: 
%   Gammaitoni, L., et al. (1998). "Stochastic resonance." Reviews of Modern Physics.
%
% 输入参数:
%   a       - 势能参数 a (标量)
%   b       - 势能参数 b (标量)
%   f       - 输入信号频率 (Hz)
%   A       - 输入信号幅值
%   D_list  - 噪声强度向量 (Noise Intensity)
%
% 输出参数:
%   snr_theory_freq - 考虑频率修正后的理论 SNR (dB)
%   snr_adiabatic   - 绝热近似理论 SNR (dB)
%   r_k             - Kramers 跃迁率

% 1. 计算势能基本参数
xm = sqrt(a / b);           % 势阱位置
dU = (a^2) / (4 * b);       % 势垒高度
omega_0 = sqrt(2 * a);      % 阱底角频率
omega_b = sqrt(a);          % 垒顶角频率

% 初始化
num_d = length(D_list);
r_k = zeros(num_d, 1);
snr_adiabatic_linear = zeros(num_d, 1);
snr_theory_linear = zeros(num_d, 1);

% 2. 迭代计算不同噪声强度下的理论值
for i = 1:num_d
    current_d = D_list(i);
    
    % Kramers 速率计算
    r_k(i) = (omega_0 * omega_b / (2 * pi)) * exp(-dU / current_d);

    % r_k(i) = (w0 * wb / (2*pi)) * exp(-dU / D_fixed);
    % snr_adiabatic(i) = (xm^2 * A^2 * r_k(i)) / (4 * D_fixed^2);
    snr_adiabatic(i) = (pi/2) * (A*xm/current_d)^2 * r_k(i); % A. 绝热近似理论
    
    % 绝热近似公式 (线性尺度)
    % SNR = (pi * A^2 * xm^2) / (2 * D^2) * r_k
    % snr_adiabatic_linear(i) = (pi * A^2 * xm^2 * r_k(i)) / (2 * current_d^2);
    
    % 频率修正因子 (Lorentzian Roll-off)
    correction_factor = 1 / (1 + (pi * f / r_k(i))^2);
    snr_theory_linear(i) = snr_adiabatic_linear(i) * correction_factor;
end

% 3. 转换为分贝 (dB)
% 注意：理论公式通常给出的是功率比，转换为dB需 10*log10
snr_theory_freq = 10 * log10(snr_theory_linear + eps);
% snr_adiabatic = 10 * log10(snr_adiabatic_linear + eps);

end