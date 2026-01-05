function [ami_score, results] = AMI2(x, fs, max_lag_sec)
% AMI2 基于驻留时间直方图的改进 AMI 指标
%
% 输入:
%   x           : 输入信号向量
%   fs          : 采样率 (Hz)
%   max_lag_sec : (可选) 最大搜索滞后时间(秒)，默认搜索到 0.5*信号长度
%
% 输出:
%   ami_score   : 最终的 SR 衡量指标 (标量)，越高代表共振越强
%   results     : 包含中间结果的结构体 (便于调试和绘图)
%       .lags       - 计算 AMI 的滞后时间轴 (秒)
%       .ami_curve  - 对应的 AMI 曲线
%       .tau_mode   - 统计得到的“主驻留时间” (秒)
%       .tau_samples- 主驻留时间对应的采样点数

% --- 1. 参数初始化与预处理 ---
x = x(:);
N = length(x);
if nargin < 3 || isempty(max_lag_sec)
    max_lag_sec = (N / 2) / fs;
end

% --- 2. 施密特触发器式二值化 (或简单二值化) ---
s = int8(x >= 0);

% --- 3. 提取驻留时间 (Pulse Widths) ---
% 找到状态翻转的索引
edges = find(diff(s) ~= 0);

if length(edges) < 2
    % 几乎没有翻转，说明信号甚至没跨过势垒，SR 未发生
    ami_score = 0;
    results.lags = [];
    results.ami_curve = [];
    results.tau_mode = 0;
    results.tau_samples = 0;
    return;
end

% 计算相邻翻转点之间的距离 (即脉冲宽度/驻留时间)
pulse_widths = [edges(1); diff(edges)]; % 包括首
% pulse_widths = diff(edges); % 单位: 采样点数

% 过滤掉极短的毛刺 (例如 < 3 点)，这些通常是高频噪声
valid_widths = pulse_widths(pulse_widths >= 3);
%
if isempty(valid_widths)
    tau_samples = round(trimean(pulse_widths, 10)); % 备选：如果全是毛刺，取平均
    fprintf('警告：所有驻留时间均过短，可能为噪声主导，取均值估计驻留时间 (%d 点)。\n', tau_samples);
else
    tau_samples = round(trimmean(valid_widths, 10)); % 直接取均值作为稳健估计
end

% 转换为物理时间 (秒)
tau_mode_sec = tau_samples / fs;
% fprintf('Estimated mode residence time: %.3f s (%d samples)\n', tau_mode_sec, tau_samples);

% --- 5. 计算 AMI 曲线 ---
% 策略：我们不计算所有 lag，重点关注 tau_mode 附近的 lag
% SR 的特征是：在 T_signal/2 (即 tau_mode) 处有很强的相关性
% 搜索范围：从 0.5 * tau_mode 到 1.5 * tau_mode
% 如果 tau_samples 很小(噪声主导)，则计算一段常规范围
if tau_samples < fs
    % 驻留时间太短，可能是纯噪声，按常规计算一小段
    ami_score = 0; results = [];
    fprintf('警告：驻留时间过短 (%d 点)，可能为噪声主导，计算常规 AMI 曲线。\n', tau_samples);
    return;
else
    % 在主峰附近搜索最大互信息
    center = tau_samples;
    width = round(tau_samples * 0.5);
    start_lag = max(1, center - width);
    end_lag = min(floor(N/2), center + width);
    search_lags = start_lag:end_lag;
end

% 确保不超过用户设定的最大 lag
max_lag_samples = round(max_lag_sec * fs);
search_lags = search_lags(search_lags <= max_lag_samples);

if isempty(search_lags)
    ami_score = 0; results = []; return;
end

ami_curve = zeros(length(search_lags), 1);

% 计算 AMI 主循环
% p0, p1 是全局概率，不需要在循环里重算
p1_global = sum(s == 1) / N;
p0_global = 1 - p1_global;
% H_s = - (p0_global * log2(p0_global + eps) + p1_global * log2(p1_global + eps));

for i = 1:length(search_lags)
    lag = search_lags(i);
    
    % 截取序列
    s1 = s(1 : N-lag);
    s2 = s(1+lag : N);
    len = length(s1);
    
    % 快速计算联合概率 (00, 01, 10, 11)
    % 利用数值技巧：s1*2 + s2 得到 0, 1, 2, 3 四种状态
    combined = s1 * 2 + s2;
    
    n00 = sum(combined == 0);
    n01 = sum(combined == 1);
    n10 = sum(combined == 2);
    n11 = sum(combined == 3);
    
    probs = [n00, n01, n10, n11] ./ len;
    
    % 互信息 I(X;Y) = sum p(x,y) * log2( p(x,y) / (p(x)*p(y)) )
    % 这里做个近似：假设局部边缘概率 p(x), p(y) 接近全局 p0, p1
    % 这样计算速度快且对于稳态信号误差极小
    mi_val = 0;
    if probs(1) > 0, mi_val = mi_val + probs(1) * log2(probs(1)/(p0_global*p0_global)); end
    if probs(2) > 0, mi_val = mi_val + probs(2) * log2(probs(2)/(p0_global*p1_global)); end
    if probs(3) > 0, mi_val = mi_val + probs(3) * log2(probs(3)/(p1_global*p0_global)); end
    if probs(4) > 0, mi_val = mi_val + probs(4) * log2(probs(4)/(p1_global*p1_global)); end
    
    ami_curve(i) = mi_val;
end

% --- 6. 提取最终指标 ---
% [peak_heights, ~] = findpeaks(ami_curve, ...
%     'MinPeakHeight', 0.01);
[peak_heights, ~] = findpeaks(ami_curve);

if isempty(peak_heights)
    ami_score = max(ami_curve) - prctile(ami_curve, 20);
else
    ami_score = max(peak_heights);
end

% --- 7. 打包结果 ---
results.lags = search_lags / fs;
results.ami_curve = ami_curve;
results.tau_mode = tau_mode_sec;
results.tau_samples = tau_samples;

end