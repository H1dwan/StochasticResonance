function fitness = SNREvaluator(theta_vec, clean_sig, noise_seq, fs, f0)
%SNREvaluator  计算 HSUBSR 系统参数下输出信号的信噪比并返回负值作为适应度
%   fitness = SNREvaluator(theta_vec, clean_sig, noise_seq, fs, f0)
%
%   输入参数：
%     theta_vec : [xm, dU, shape_factor] 参数向量，依次用于校准 HSUBSR 势函数系数
%     clean_sig : 原始无噪输入信号向量
%     noise_seq : 与输入同长度的噪声序列
%     fs        : 采样频率 (Hz)
%     f0        : 输入信号的目标频率 (Hz)，用于 SNR 计算
%
%   输出参数：
%     fitness   : 适应度值，取输出 SNR 的相反数以便于最小化型优化算法
%

% 0. 解包待优化参数
[a, b, k1, k2] = CalibrateHSUBSR(theta_vec(1), theta_vec(2), theta_vec(3));

% 1. 仅截取稳态段用于指标计算，剔除前 10% 过渡期
n_samples    = length(clean_sig);
steady_start = round(0.1 * n_samples);

% 2. 构建漂移函数并数值积分得到输出轨迹
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x_out = RK4Solver(drift_func, clean_sig + noise_seq, 1/fs);
x_out = x_out - mean(x_out); % 去直流分量
x_steady = x_out(steady_start+1:end);

% 3. 计算输出信号的 SNR 作为适应度指标
snr_val = SNRo2(x_steady, fs, f0);
if isinf(snr_val)
    snr_val = -100; % 对于无穷大的 SNR，设置为 -100 以避免优化问题
end

fitness = -snr_val;

end