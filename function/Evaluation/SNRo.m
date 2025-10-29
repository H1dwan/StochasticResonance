function snr = SNRo(s_out, fs, f_target)
% SNRo 计算随机共振系统输出信号的信噪比(SNR)
% 
% 该函数通过FFT分析输出信号，在频域中计算目标频率处的信号功率，
% 并根据公式计算输出信噪比。
%输入参数:
%   s_out     - 输出信号向量
%   fs        - 采样频率(Hz)
%   f_target  - 目标信号频率(Hz)
%
% 返回值:
%   snr       - 输出信噪比(dB)
%
% 参考公式:
%   公式(2-54): NN = sum of squared values of first half of signal
%   公式(2-55): SNR = 10*log10(SP/(NN-SP))

    % 输入参数验证 - 修复：添加参数有效性检查
    if ~isnumeric(s_out) || isempty(s_out)
        error('SNRo:InvalidSignal', '输入信号必须是非空数值向量');
    end
    
    if ~isnumeric(fs) || fs <= 0
        error('SNRo:InvalidSamplingRate', '采样频率必须是正数');
    end
    
    if ~isnumeric(f_target) || f_target <= 0 || f_target > fs/2
        error('SNRo:InvalidTargetFrequency', '目标频率必须在(0, fs/2]范围内');
    end

    % 计算输出信号的频谱特性
    N = length(s_out);
    % 修复：去除均值以减少直流分量影响
    s_out = s_out - mean(s_out);
    % 修复：使用窗函数减少频谱泄漏
    window = hann(N);
    s_out = s_out .* window;
    
    % 修复：使用2的幂次长度进行FFT以提高效率和频率分辨率
    N_fft = 2^nextpow2(N);
    X = fft(s_out, N_fft);
    X_mag = abs(X);
    % 修复：频谱幅度归一化
    X_mag = X_mag / N;

    % 提取正频率部分以减少计算量
    half_N = floor(N_fft/2) + 1;
    X_mag_half = X_mag(1:half_N);

    % 构建对应的频率轴
    frequencies = (0:half_N-1) * fs / N_fft;

    % 定位目标频率在频谱中的位置
    [~, idx] = min(abs(frequencies - f_target));
    
    SP = sum(X_mag_half(idx).^2);
    
    NN = sum(X_mag_half.^2);

    % 根据公式计算输出信噪比
    if NN <= 0
        warning('SNRo:InvalidNoisePower', '噪声功率估计过小或为负');
        snr = Inf;
    else
        snr = 10 * log10(SP / (NN-SP));
    end
end