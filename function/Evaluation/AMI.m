function [ami, valid_curve] = AMI(x, fs, k)
% AMI:
% 输入: x - 信号; fs - 采样率
% 输出: AMI_pure - 纯互信息波动

x = x(:);
N = length(x);

% 理论依据：SR中粒子在势阱间的跳变是大尺度的，阱内抖动是小尺度的。
% 利用信号的统计分布（标准差）来区分这两种尺度。
sigma = std(x);
mu = mean(x);

% 设定滞后阈值：均值上下 0.3 倍标准差 (经验值，适应性强)
th_high = mu + 0.3 * sigma;
th_low  = mu - 0.3 * sigma;

s = zeros(N, 1);
state = 1; % 初始状态
if x(1) < mu, state = 0; end

for i = 1:N
    if x(i) > th_high
        state = 1;
    elseif x(i) < th_low
        state = 0;
    end
    s(i) = state;
end
% 此时 s 是去除了毛刺的纯净二值序列

% --- 1. AMI 计算核心参数 ---
max_lag = min(round(k * fs), floor(N/2)); % 最大延迟取 k 秒
% s = int8(x >= 0); % 状态符号化 (0 或 1)

AMI_curve = zeros(max_lag, 1);

% 预计算全局概率
p1 = sum(s) / N;
p0 = 1 - p1;
if p1 < 1e-3 || p0 < 1e-3
    ami = 0; valid_curve = zeros(max_lag, 1); return;
end

% --- 2. 循环计算延迟互信息曲线 ---
for tau = 1:max_lag
    s_now = s(1 : N-tau);
    s_fut = s(1+tau : N);
    len = length(s_now);
    
    % 联合概率统计
    n11 = sum(s_now & s_fut); n00 = sum(~s_now & ~s_fut);
    n01 = sum(~s_now & s_fut); n10 = len - n11 - n00 - n01;
    
    probs = [n00, n01, n10, n11] / len;
    
    % 边际概率
    pi_0 = (n00 + n01) / len; pi_1 = 1 - pi_0;
    pj_0 = (n00 + n10) / len; pj_1 = 1 - pj_0;
    
    marg_i = [pi_0, pi_0, pi_1, pi_1];
    marg_j = [pj_0, pj_1, pj_0, pj_1];
    
    % 计算互信息
    mi_val = 0;
    if probs(1)>eps, mi_val = mi_val + probs(1)*log2(probs(1)/(marg_i(1)*marg_j(1))); end
    if probs(2)>eps, mi_val = mi_val + probs(2)*log2(probs(2)/(marg_i(2)*marg_j(2))); end
    if probs(3)>eps, mi_val = mi_val + probs(3)*log2(probs(3)/(marg_i(3)*marg_j(3))); end
    if probs(4)>eps, mi_val = mi_val + probs(4)*log2(probs(4)/(marg_i(4)*marg_j(4))); end
    
    AMI_curve(tau) = mi_val;
end

cutoff_idx = max(1, floor(length(AMI_curve) * 0.2));
valid_curve = AMI_curve(cutoff_idx:end);

if isempty(valid_curve)
    ami = 0;
else
    % 纯AMI指标 = 曲线的标准差
    % ami = std(valid_curve);
    baseline = prctile(valid_curve, 20);
    peak_val = max(valid_curve);
    ami = peak_val - baseline;
    % ami = mean(valid_curve);
end
end