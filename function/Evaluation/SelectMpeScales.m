function [scales, info] = SelectMpeScales(x, fs, m)
% 盲选MPE尺度：基于两态化驻留时间 + 统计长度约束
%
% 输入：
%   x   : 经过 SR 系统的输出信号（列向量）
%   fs  : 采样频率（Hz）
%   m   : MPE嵌入维度
%
% 输出：
%   scales : 选定的尺度列表
%   info   : 选尺度的辅助信息结构体

x = x(:);
n = length(x);

% ---- 1) 估计 |x| 的众数（粗略用直方图峰值代替）----
abs_x = abs(x);
nbins = 100;
[counts, edges] = histcounts(abs_x, nbins);
[~, idx] = max(counts);
mode_abs = 0.5 * (edges(idx) + edges(idx+1));

% ---- 2) 施密特触发两态化，减少抖动翻转 ----
theta = 0.7 * mode_abs;       % 推荐：0.6~0.8 * mode(|x|)
state = 0;
states = zeros(n,1);

for i = 1:n
    xi = x(i);
    if state == 0
        if xi >= theta
            state = 1;
        elseif xi <= -theta
            state = -1;
        end
    elseif state == 1
        if xi <= -theta
            state = -1;
        end
    else % state == -1
        if xi >= theta
            state = 1;
        end
    end
    states(i) = state;
end

% 将0状态用前值填充
for i = 2:n
    if states(i) == 0
        states(i) = states(i-1);
    end
end
% if states(1) == 0
%     states(1) = 1;
% end
% 首段全为0时，用首个非零状态回填前缀
if any(states == 0)
    firstNz = find(states ~= 0, 1);
    if ~isempty(firstNz)
        states(1:firstNz) = states(firstNz);
    end
end

% ---- 3) 计算驻留时间分布 ----
% dur_samples = [];
% cur = states(1);
% len = 1;
% for i = 2:n
%     if states(i) == cur
%         len = len + 1;
%     else
%         dur_samples(end+1,1) = len; %#ok<AGROW>
%         cur = states(i);
%         len = 1;
%     end
% end
% dur_samples(end+1,1) = len;
edges = find(diff(states) ~= 0);
if isempty(edges)
    % 没有翻转，说明信号没跨过势垒
    dur_samples = n;
else
    dur_samples = [edges(1); diff(edges)];
end
dur_sec = dur_samples / fs;

t_res_mean   = mean(dur_sec);
t_res_median = median(dur_sec);
t_res_trim   = trimmean(dur_sec, 10); % 去掉10%极值后的均值

% ---- 4) 由驻留时间给物理上限 + 由长度给统计上限 ----
s_min = max(1, floor(t_res_trim * fs / (6*m-1)));   % 绝热低频下建议至少从1秒尺度起步(10点)
s_max_phy  = max(1, floor(t_res_trim * fs / (10)));

C = 6; % 推荐：m=3时取6（平均每种序模式≈6次以上）
s_max_stat = floor(n / (C * factorial(m)));

s_max = min(s_max_phy, s_max_stat);
if s_min > s_max
    s_min = 1;
    s_max = 2;
end
scales = (s_min:4:s_max).';
if length(scales) <= 2
    scales = (s_min:2:s_max).';
end

if s_min == 1
    scales = 1;
end

info = struct();
info.mode_abs = mode_abs;
info.theta = theta;
info.t_res_mean = t_res_mean;
info.t_res_median = t_res_median;
info.t_res_trim = t_res_trim;
info.s_min = s_min;
info.s_max_phy = s_max_phy;
info.s_max_stat = s_max_stat;
info.s_max = s_max;

end