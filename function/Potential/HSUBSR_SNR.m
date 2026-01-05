function [snr_theory_freq, snr_adiabatic, r_k] = HSUBSR_SNR(xm, dU, shape, f, A, D_list)


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
    snr_adiabatic(i) = (xm^2 * A^2 * r_k(i)) / (4 * D_fixed^2);
    correction = 1 / (1 + (pi * f / r_k(i))^2); % 频率修正因子
    snr_theory_freq(i) = snr_adiabatic(i) * correction;
end

end