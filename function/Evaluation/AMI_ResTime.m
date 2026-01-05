function [ami, valid_curve, tau_valid] = AMI_ResTime(x, fs, k, tau_factor)
% AMI_ResTime  基于驻留时间自适应裁剪小延迟的改进 AMI 指标
%
%   [ami, valid_curve, tau_valid] = AMI_ResTime(x, fs, k)
%   [ami, valid_curve, tau_valid] = AMI_ResTime(x, fs, k, tau_factor)
%
% 输入参数：
%   x          : 输入信号（一维向量，行/列均可）
%   fs         : 采样频率 (Hz)
%   k          : 最大互信息滞后时间窗口 (秒)，即 tau_max = k
%   tau_factor : (可选) 驻留时间缩放因子，默认 0.5
%                tau_min = tau_factor * tau_res
%
% 输出参数：
%   ami        : 改进后的 AMI 指标（标量）
%   valid_curve: 剪裁后的 AMI(tau) 曲线（从 tau_min 起）
%   tau_valid  : valid_curve 对应的滞后时间轴 (秒)，便于画图
%
% 理论说明：
%   - SR 中我们关注的是“跨阱跃迁随时间的相关性”，其时间尺度由驻留时间 tau_res 决定；
%   - 对 tau << tau_res 的延迟，互信息主要反映的是“短时自相关/记忆”，会导致在小噪声 D 时
%     AMI 曲线偏高，从而把最优 D 向左拉偏；
%   - 本函数先用 ZCR(x, fs) 估计平均驻留时间 tau_res ≈ 1/ZCR，再令
%       tau_min = tau_factor * tau_res
%     只在 tau ≥ tau_min 的延迟上统计 AMI，尽可能抑制短时自相关的影响。
%

% -------------------- 0. 输入整理与默认参数 ------------------------
if nargin < 4 || isempty(tau_factor)
    tau_factor = 0.5;   % 默认只保留大约 ≥ 0.5 * tau_res 的延迟
end

x = x(:);                          % 转列向量
num_samples = length(x);

if num_samples < 10
    % 样本太少，直接返回 0
    ami = 0;
    valid_curve = [];
    if nargout >= 3
        tau_valid = [];
    end
    return;
end

% -------------------- 1. 滞回二值化（抑制阱内抖动） ----------------
s = int8(x >= 0); % 状态符号化 (0 或 1)

% -------------------- 2. 利用 ZCR 估计驻留时间 tau_res -------------
%
% ZCR(x, fs) 返回“每秒翻转次数”，近似可视为跃迁率 r_flip，
% 则平均驻留时间 tau_res ≈ 1 / r_flip。
try
    zcr_val = ZCR(x, fs);   % 若函数在路径中，则直接使用
    % zcr_val = EstimateEscapeRate(x, fs, 1);
catch
    warning('AMI_ResTime:ZCRNotFound', ...
        '未找到 ZCR.m，使用退化估计：假定整个序列只含一个驻留段。');
    zcr_val = 0;
end

if zcr_val <= 0
    % 极端情况：几乎没有跨阈值翻转，用整段观测时间代表驻留时间
    tau_res_sec = num_samples / fs;
else
    tau_res_sec = 1.0 / zcr_val;
end

% fprintf('Estimated residence time tau_res = %.3f s (ZCR = %.3f Hz)\n', ...
%     tau_res_sec, zcr_val);

% -------------------- 3. 定义最大滞后与 tau_min -------------------
max_lag = min(round(k * fs), floor(num_samples / 2));
if max_lag < 1
    ami = 0;
    valid_curve = [];
    if nargout >= 3
        tau_valid = [];
    end
    return;
end

tau_vec = (1:max_lag).' / fs;      % 所有滞后点对应的时间 (秒)

% 自适应小延迟下限：tau_min = tau_factor * tau_res
tau_min_sec = tau_factor * tau_res_sec;
idx_start = find(tau_vec >= tau_min_sec, 1, 'first');

% 若 tau_min 过大导致无有效点，则退回到经验 20% 截断
if isempty(idx_start)
    idx_start = max(1, floor(0.2 * max_lag));
end

% fprintf('Using tau_min = %.3f s (index %d of %d)\n', ...
%     tau_vec(idx_start), idx_start, max_lag);

% -------------------- 4. 计算二值序列的延迟互信息 I(tau) ---------
ami_curve = zeros(max_lag, 1);

for lag = 1:max_lag
    s1 = s(1:(num_samples - lag));       % 当前时刻
    s2 = s((1 + lag):num_samples);       % 滞后 lag 后的时刻
    pair_len = length(s1);
    
    if pair_len <= 0
        ami_curve(lag) = 0;
        continue;
    end
    
    % 统计 2x2 联合频数 N_ij
    n00 = sum((s1 == 0) & (s2 == 0));
    n01 = sum((s1 == 0) & (s2 == 1));
    n10 = sum((s1 == 1) & (s2 == 0));
    n11 = sum((s1 == 1) & (s2 == 1));
    
    p00 = n00 / pair_len;
    p01 = n01 / pair_len;
    p10 = n10 / pair_len;
    p11 = n11 / pair_len;
    
    % 边缘概率：p_i = P(s1 = i), q_j = P(s2 = j)
    p0 = p00 + p01;
    p1 = p10 + p11;
    q0 = p00 + p10;
    q1 = p01 + p11;
    
    % 互信息 I = sum_ij p_ij log2( p_ij / (p_i * q_j) )
    mi_val = 0;
    
    if p00 > 0 && p0 > 0 && q0 > 0
        mi_val = mi_val + p00 * log(p00 / (p0 * q0));
    end
    if p01 > 0 && p0 > 0 && q1 > 0
        mi_val = mi_val + p01 * log(p01 / (p0 * q1));
    end
    if p10 > 0 && p1 > 0 && q0 > 0
        mi_val = mi_val + p10 * log(p10 / (p1 * q0));
    end
    if p11 > 0 && p1 > 0 && q1 > 0
        mi_val = mi_val + p11 * log(p11 / (p1 * q1));
    end
    
    % 换算为 bit（log 默认是自然对数）
    ami_curve(lag) = mi_val / log(2);
end

% -------------------- 5. 剪裁掉 tau < tau_min 的小延迟 ------------
valid_curve = ami_curve(idx_start:end);
if nargout >= 3
    tau_valid = tau_vec(idx_start:end);
end

% -------------------- 6. 标量 AMI 指标：峰值-基线 ------------------
if isempty(valid_curve)
    ami = 0;
else
    [pks, locs] = findpeaks(valid_curve, 'MinPeakProminence', 0.01);
    if isempty(pks)
        ami = 0; % 无相关性
    else
        ami = pks(1);
    end
    % % ---- 2) 定义“接近 0”的阈值（绝对 + 相对） ----
    % peak0      = max(valid_curve);
    % thr_abs    = 0.01;                % 绝对阈值：0.001 bit，可调
    % thr_rel    = 0.05 * peak0;         % 相对阈值：峰值的 1%
    % thr_zero   = max(thr_abs, thr_rel);  % 取更大的那一个，保守一些
    
    % % ---- 3) 要求有一段连续点都低于 thr_zero，表示进入“近 0 平台” ----
    % below      = valid_curve < thr_zero;
    % run_len    = max(3, round(0.02 * length(valid_curve)));  % 至少 2% 长度的连续区
    % below_int  = conv(double(below), ones(run_len,1), 'same');
    
    % idx_split  = find(below_int >= run_len, 1, 'first');
    
    % if isempty(idx_split)
    %     % 找不到“近 0 平台”，说明相关性衰减很慢或数据太短
    %     % 可以退回一个经验分割点，比如 60% 位置
    %     idx_split = max(1, floor(0.6 * length(valid_curve)));
    % end
    
    % % ---- 4) 在“近 0 后的振荡区”上再计算指标 ----
    % search_curve = valid_curve(idx_split:end);
    % if isempty(search_curve)
    %     ami = 0;
    % else
    %     % 基线用该区间的 20% 分位数
    %     baseline = prctile(search_curve, 20);
    %     peak_val = max(search_curve);
    %     ami = peak_val - baseline;
    % end
end
end