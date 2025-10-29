function snr = calculate_SNR(signal, fs, f_target, varargin)
% calculate_SNR 计算信号在特定频率处的信噪比(SNR)
% 
% 该函数计算给定信号在指定特征频率处的信噪比，可以通过多种方法估算噪声
% 水平，并提供信号预处理选项。
% 
% 输入参数:
%   signal    - 输入信号向量
%   fs        - 采样频率(Hz)
%   f_target  - 目标信号频率(Hz)
%   varargin  - 可选参数（键值对形式）:
%     'NoiseEstMethod' - 噪声估算方法：'spectrum_subtraction' (默认), 
%                        'bandwidth_estimation', 'reference_noise'
%     'SignalBandwidth' - 信号带宽(Hz)，默认为目标频率的5%
%     'ReferenceNoise'  - 参考噪声信号（当NoiseEstMethod为'reference_noise'时使用）
%     'Preprocess'      - 是否对信号进行预处理，默认为true
%     'Window'          - 窗函数类型，默认为'hann'
%     'OutputUnit'      - 输出单位：'dB'(默认)或'linear'
% 
% 返回值:
%   snr       - 信噪比（dB或线性比例）
% 
% 示例:
%   % 基本用法
%   snr = calculate_SNR(signal, 1000, 50);
%   
%   % 使用带宽估计法计算SNR
%   snr = calculate_SNR(signal, 1000, 50, 'NoiseEstMethod', 'bandwidth_estimation', ...
%                       'SignalBandwidth', 2);
%   
%   % 使用参考噪声计算SNR
%   snr = calculate_SNR(signal, 1000, 50, 'NoiseEstMethod', 'reference_noise', ...
%                       'ReferenceNoise', noise_signal);

    % 解析默认参数
    p = inputParser;
    addRequired(p, 'signal', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'fs', @(x) isnumeric(x) && x > 0);
    addRequired(p, 'f_target', @(x) isnumeric(x) && x > 0 && x <= fs/2);
    addParameter(p, 'NoiseEstMethod', 'spectrum_subtraction', ...
                 @(x) ismember(lower(x), {'spectrum_subtraction', 'bandwidth_estimation', 'reference_noise'}));
    addParameter(p, 'SignalBandwidth', f_target * 0.05, ...
                 @(x) isnumeric(x) && x > 0);
    addParameter(p, 'ReferenceNoise', [], @isnumeric);
    addParameter(p, 'Preprocess', true, @islogical);
    addParameter(p, 'Window', 'hann', @ischar);
    addParameter(p, 'OutputUnit', 'dB', ...
                 @(x) ismember(lower(x), {'db', 'linear'}));
    
    % 解析输入参数
    parse(p, signal, fs, f_target, varargin{:});
    
    % 提取参数值
    signal = p.Results.signal;
    fs = p.Results.fs;
    f_target = p.Results.f_target;
    noise_method = lower(p.Results.NoiseEstMethod);
    signal_bw = p.Results.SignalBandwidth;
    ref_noise = p.Results.ReferenceNoise;
    preprocess = p.Results.Preprocess;
    window_type = p.Results.Window;
    output_unit = lower(p.Results.OutputUnit);
    
    % 参数验证
    if f_target > fs/2
        error('目标频率不能超过奈奎斯特频率（fs/2）');
    end
    
    if strcmp(noise_method, 'reference_noise') && isempty(ref_noise)
        error('使用参考噪声法时必须提供参考噪声信号');
    end
    
    % 信号预处理
    if preprocess
        % 去除均值
        signal = signal - mean(signal);
        
        % 应用窗函数
        window = window(@(n) eval([window_type '(n)']), length(signal));
        signal = signal .* window';
    end
    
    % 计算信号频谱
    N = length(signal);
    if N < 1024
        % 补零以提高频率分辨率
        N_fft = 2^nextpow2(1024);
    else
        N_fft = 2^nextpow2(N);
    end
    
    X = fft(signal, N_fft);
    X_mag = abs(X) / N;  % 归一化幅度
    X_mag_half = X_mag(1:N_fft/2+1);
    frequencies = (0:N_fft/2) * fs / N_fft;
    
    % 计算信号功率
    switch noise_method
        case 'spectrum_subtraction'
            % 频谱减法：从总功率中减去目标频率处的功率
            [~, f_idx] = min(abs(frequencies - f_target));
            
            % 计算信号带宽内的总功率
            bw_lower = f_target - signal_bw/2;
            bw_upper = f_target + signal_bw/2;
            
            % 找到带宽范围内的频率索引
            bw_idx = frequencies >= bw_lower & frequencies <= bw_upper;
            
            % 计算信号功率（频率域功率）
            signal_power = sum(X_mag_half(bw_idx).^2);
            
            % 计算总功率
            total_power = sum(X_mag_half.^2);
            
            % 估算噪声功率
            noise_power = total_power - signal_power;
            
        case 'bandwidth_estimation'
            % 带宽估计法：只考虑目标频率附近的信号，其余视为噪声
            [~, f_idx] = min(abs(frequencies - f_target));
            
            % 计算信号带宽内的总功率
            bw_lower = f_target - signal_bw/2;
            bw_upper = f_target + signal_bw/2;
            
            % 找到带宽范围内的频率索引
            signal_idx = frequencies >= bw_lower & frequencies <= bw_upper;
            
            % 计算信号功率
            signal_power = sum(X_mag_half(signal_idx).^2);
            
            % 计算噪声功率（非信号带宽内的功率）
            noise_idx = ~signal_idx;
            noise_power = sum(X_mag_half(noise_idx).^2);
            
        case 'reference_noise'
            % 参考噪声法：使用参考噪声信号估算噪声水平
            if preprocess
                % 对参考噪声也进行相同的预处理
                ref_noise = ref_noise - mean(ref_noise);
                window = window(@(n) eval([window_type '(n)']), length(ref_noise));
                ref_noise = ref_noise .* window';
            end
            
            % 计算参考噪声的频谱
            N_noise = length(ref_noise);
            if N_noise < N_fft
                X_noise = fft(ref_noise, N_fft);
            else
                X_noise = fft(ref_noise(1:N_fft));
            end
            
            X_noise_mag = abs(X_noise) / N_noise;
            X_noise_mag_half = X_noise_mag(1:N_fft/2+1);
            
            % 计算信号在目标频率处的功率
            [~, f_idx] = min(abs(frequencies - f_target));
            
            % 计算信号带宽内的总功率
            bw_lower = f_target - signal_bw/2;
            bw_upper = f_target + signal_bw/2;
            
            % 找到带宽范围内的频率索引
            bw_idx = frequencies >= bw_lower & frequencies <= bw_upper;
            
            % 计算信号功率
            signal_power = sum(X_mag_half(bw_idx).^2);
            
            % 计算相同带宽内的噪声功率
            noise_power = sum(X_noise_mag_half(bw_idx).^2);
    end
    
    % 计算信噪比
    if noise_power <= 0
        warning('噪声功率估算值过小或为负，可能导致结果不准确');
        snr = Inf;
    else
        if strcmp(output_unit, 'db')
            snr = 10 * log10(signal_power / noise_power);
        else
            snr = signal_power / noise_power;
        end
    end
end