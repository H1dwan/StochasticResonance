clc; clear; close all;

%% 1. 仿真参数设置
fs  = 5;              % 采样率
T   = 2000;           % 信号时长
N   = fs * T;
t   = (0:N-1)'/fs;

% 弱信号参数
A0  = 0.1;
f0  = 0.01;
clean_sig = A0 * sin(2*pi*f0*t);

% SR 系统参数 (a=1, b=1)
[a, b, k1, k2] = CalibrateHSUBSR(1, 0.25, 1.01);
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
% drift_func = @(x) CBSR_Dynamics(x, 1, 1);

% 扫描范围
D_list  = 0.05:0.01:0.45;   % 噪声强度范围
% D_list = 0.12;
n_D     = length(D_list);
num_mc  = 100;              % 蒙特卡洛次数，可按需调整

mpe_list2 = zeros(n_D, 1); % 存储每次MC的MPE结果
mpe_list4 = zeros(n_D, 1);
mpe_list6 = zeros(n_D, 1);

snr_list = zeros(n_D, 1);

h = 1 /fs;
%% 2. 循环计算
m = 3;  % MPE嵌入维度
parfor i = 1 : n_D
    D = D_list(i);
    mpe_rep2 = zeros(num_mc, 1);
    mpe_rep4 = zeros(num_mc, 1);
    mpe_rep6 = zeros(num_mc, 1);
    snr_rep = zeros(num_mc, 1);

    for mc_iter = 1:num_mc
        % 1. 生成含噪信号
        noise = sqrt(2*D*fs) * randn(N, 1);
        input = clean_sig + noise;
        
        % 2. RK4 求解 SR 输出
        x = RK4Solver2(drift_func, clean_sig, noise, fs);
        
        % 3. 预处理 (去瞬态 + 归一化)
        x_steady = x(round(0.1*N)+1:end);

        out2 = SR_AcfAdaptiveScale(x_steady, fs, 'm', m, 'C', 50);
        [~, mpe2, ~, ~]  = MultiScalePermEn(x_steady, out2.S, 'm', m);
        H_min2 = mean(mpe2);

        % out4 = SR_AcfAdaptiveScale(x_steady, fs, 'm', m, 'C', 4);
        % [~, mpe4, ~, ~]  = MultiScalePermEn(x_steady, out4.S, 'm', m);
        % H_min4 = mean(mpe4);

        % out6 = SR_AcfAdaptiveScale(x_steady, fs, 'm', m, 'C', 6);
        % [~, mpe6, ~, ~]  = MultiScalePermEn(x_steady, out6.S, 'm', m);
        % H_min6 = mean(mpe6);
        
        snr_rep(mc_iter) = SNRo2(x_steady, fs, f0);

        % 累积求均值/方差
        mpe_rep2(mc_iter) = H_min2;
        % mpe_rep4(mc_iter) = H_min4;
        % mpe_rep6(mc_iter) = H_min6;
    end

    mpe_list2(i) = mean(mpe_rep2);
    % mpe_list4(i) = mean(mpe_rep4);
    % mpe_list6(i) = mean(mpe_rep6);
    snr_list(i) = mean(snr_rep);
    
    if mod(i, 5) == 0
        fprintf('  D index %d/%d 完成\n', i, n_D);
    end
end

%% 绘制
figure;
yyaxis left;
plot(D_list, smooth(snr_list, 1), '-o', 'LineWidth', 2);
yyaxis right;
plot(D_list, smooth(mpe_list2, 1), '-^', 'LineWidth', 2); hold on;
% plot(D_list, smooth(mpe_list4, 1), '-s', 'LineWidth', 2);
% plot(D_list, smooth(mpe_list6, 1), '-d', 'LineWidth', 2);

results.input.fs = fs;
results.input.T = T;
results.input.A0 = A0;
results.input.f0 = f0;
results.input.clean_sig = clean_sig;
results.potential.xm = 1;
results.potential.dU = 0.25;
results.potential.shape = 1.01;
results.output.snr_list = snr_list;
results.output.mpe_list2 = mpe_list2;
results.output.mpe_list4 = mpe_list4;
results.output.mpe_list6 = mpe_list6;

%% 辅助函数
function out = SR_AcfAdaptiveScale(y, fs, varargin)
%SR_AcfAdaptiveScale  基于输出自相关(ACF)的自适应主时间尺度估计与候选尺度集合构造
%
% 【用途】
%   针对随机共振(SR)输出序列 y，仅使用输出的时间域统计结构（ACF）
%   自适应估计主时间尺度 k*（以“样本点滞后”为单位），并将 k* 映射为
%   多尺度粗粒化的候选尺度集合 S（用于后续MPE/其他指标的尺度选择）。
%
%   注意：本函数“不计算排列熵/不进行尺度优选”，仅输出：
%         - 主滞后 kStar（主时间尺度）
%         - 候选尺度集合 S
%         - 以及调试信息（ACF曲线、峰值、阈值等）
%
% 【核心思想】
%   1) 计算归一化自相关 ACF：r[k] = corr(y[n], y[n+k])
%   2) 在正滞后范围内寻找“第一个显著正的局部最大峰”作为主滞后 k*
%      - 显著性阈值：theta = 1.96/sqrt(N)，用于过滤噪声假峰
%      - 若无显著峰：退化为全局最大峰对应的滞后
%   3) 映射到尺度：s0 = floor(k*/alpha)，并在 [floor(s0/2), floor(2*s0)] 形成候选集合
%      同时施加 N/s >= Nmin 的约束，保证后续多尺度统计的有效长度
%
% 【输入】
%   y  : SR输出序列（向量）
%   fs : 采样率(Hz)。本函数内部不强制使用 fs（保留用于解释与记录）。
%
% 【可选参数】（Name-Value）
%   'kmaxCap'        : ACF最大搜索滞后上限，默认 1000
%   'kmin'           : 忽略很小滞后（抑制惯性/局部平滑引起的小峰），默认 2
%   'alpha'          : k* -> s0 映射系数，默认 4（理解为“每周期保留约alpha个粗粒化点”）
%   'Nmin'           : 粗粒化后最小有效长度约束，默认 500（要求 N/s >= Nmin）
%   'removeTransient': 丢弃前段比例（0~0.5），默认 0（不丢弃）
%   'detrendMode'    : 去趋势方式：'none'|'linear'|'diff'，默认 'none'
%   'minPeakProm'    : 峰突出度阈值（可选增强鲁棒性），默认 []（不启用）
%   'minPeakDist'    : 峰最小间隔（样本点），默认 2
%
% 【输出 out 结构体】
%   out.kStar        : ACF估计主滞后（样本点）
%   out.theta        : ACF显著阈值
%   out.s0           : 中心尺度（由k*映射得到）
%   out.S            : 候选尺度集合（列向量）
%   out.debug        : 调试信息（acfPos, kmax, peaks等）
%
% 【示例】
%   out = SR_AcfAdaptiveScale(y, 5, 'alpha', 4, 'Nmin', 500);
%   fprintf('k*=%d, s0=%d, S=[%d..%d]\n', out.kStar, out.s0, out.S(1), out.S(end));
%
% -------------------------------------------------------------------------

%% ============== 参数解析 ==============
p = inputParser;
p.addRequired('y', @(x)isnumeric(x) && isvector(x) && numel(x) > 50);
p.addRequired('fs', @(x)isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('m', 3, @(x)isnumeric(x) && isscalar(x) && x>=3);
p.addParameter('C', 2, @(x)isnumeric(x) && isscalar(x) && x>=2);
p.addParameter('kmaxCap', 1000, @(x)isnumeric(x) && isscalar(x) && x>=10);
p.addParameter('kmin', 5, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('alpha', 4, @(x)isnumeric(x) && isscalar(x) && x>=2);
p.addParameter('Nmin', 500, @(x)isnumeric(x) && isscalar(x) && x>=50);

p.addParameter('removeTransient', 0, @(x)isnumeric(x) && isscalar(x) && x>=0 && x<0.5);
p.addParameter('detrendMode', 'none', @(x)ischar(x) || isstring(x));

p.addParameter('minPeakProm', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
p.addParameter('minPeakDist', 2, @(x)isnumeric(x) && isscalar(x) && x>=1);

p.parse(y, fs, varargin{:});
opt = p.Results;

y = y(:);            % 强制列向量
N0 = numel(y);

% ============== 可选：丢弃前段瞬态 ==============
if opt.removeTransient > 0
    nDrop = round(opt.removeTransient * N0);
    if nDrop >= 1 && (N0 - nDrop) > 50
        y = y(nDrop+1:end);
    end
end

% ============== 可选：去趋势处理 ==============
switch lower(string(opt.detrendMode))
    case "none"
        % 不处理
    case "linear"
        y = detrend(y);  % 线性去趋势
    case "diff"
        y = diff(y);     % 一阶差分（强去趋势，可能改变结构，慎用）
    otherwise
        error('detrendMode 仅支持: none | linear | diff');
end

% ============== 标准化（避免幅值影响ACF数值） ==============
y = y - mean(y);
sd = std(y);
if sd < eps
    error('输入序列方差过小，无法计算ACF。');
end
x = y / sd;
N = numel(x);

%% ============== 1) 计算ACF并估计主滞后 k* ==============
% 选择搜索最大滞后：
%   - kmaxCap 控制计算量
%   - N/10 避免远滞后估计不稳定
kmax = min([opt.kmaxCap, floor(N/10), N-2]);
if kmax < (opt.kmin + 2)
    % 数据太短或参数不合理，兜底输出
    out = local_default_out(fs);
    return;
end

% 归一化自相关：xcorr 返回 [-kmax..kmax]，取正滞后 0..kmax
acfAll = xcorr(x, kmax, 'coeff');
acfPos = acfAll(kmax+1:end);       % lag = 0..kmax
acfPos(1) = 0;                     % 忽略 lag=0（必为最大值）

% ACF显著阈值（经验）：白噪声非零滞后ACF约在 ±1.96/sqrt(N) 内
theta = 1.96/sqrt(N);

% 峰值搜索：可选增加 MinPeakProminence 增强抗假峰能力
if isempty(opt.minPeakProm)
    [pks, locs] = findpeaks(acfPos, 'MinPeakDistance', opt.minPeakDist);
else
    [pks, locs] = findpeaks(acfPos, ...
        'MinPeakDistance', opt.minPeakDist, ...
        'MinPeakProminence', opt.minPeakProm);
end

% 候选峰：lag >= kmin 且峰值 > theta
cand = find(locs >= (opt.kmin+1) & pks > theta);

if ~isempty(cand)
    % 规则：取“第一个显著正峰”作为主滞后 k*
    % 目的：避免倍周期峰（2T/3T）或偶然更高峰导致 k* 偏大
    loc   = locs(cand(1));
    kStar = loc - 1;     % lag = index - 1
    modeSel = "first_significant_peak";
else
    % 退化：在 [kmin..kmax] 内选全局最大ACF对应滞后
    [~, idx] = max(acfPos((opt.kmin+1):end));
    kStar = opt.kmin + idx;
    modeSel = "global_max_fallback";
end

% ============== 2) k* -> 候选尺度集合 S ==============
% 映射：s0 = floor(k*/alpha)
% 含义：若 k* 近似表示一个周期的采样点数，则希望粗粒化后每周期保留约 alpha 个点
s0 = max(1, floor(kStar / (opt.C + opt.m - 1)));

% 候选集合：围绕 s0 做局部搜索，而非全尺度盲扫
sLo = max(1, floor(s0/2));
sHi = max(sLo, floor(2*s0));

% 约束：粗粒化后有效长度至少 Nmin -> N/s >= Nmin -> s <= floor(N/Nmin)
sHi = min(sHi, floor(N / opt.Nmin));

% 兜底处理：若约束导致区间为空，给一个保守范围
if sHi < sLo
    sLo = 1;
    sHi = min(max(1, floor(N/opt.Nmin)), 20);
end

S = (sLo:sHi)';
if isempty(S)
    S = 1;
end

% ============== 输出 ==============
out = struct();
out.kStar = kStar;
out.theta = theta;
out.s0    = s0;
out.S     = S;

% 调试信息（用于论文画图或检查）
dbg = struct();
dbg.fs        = fs;
dbg.N         = N;
dbg.kmax      = kmax;
dbg.kmin      = opt.kmin;
dbg.alpha     = opt.alpha;
dbg.Nmin      = opt.Nmin;
dbg.acfPos    = acfPos(:);
dbg.pks       = pks(:);
dbg.locs      = locs(:);           % loc index (lag = loc-1)
dbg.modeSel   = modeSel;           % 选择模式（首显著峰/退化最大峰）
dbg.sLo       = sLo;
dbg.sHi       = sHi;
out.debug     = dbg;
end

function out = local_default_out(fs)
% 兜底输出（当数据太短或kmax设置不合理时）
out = struct();
out.kStar = 1;
out.theta = NaN;
out.s0    = 1;
out.S     = 1;
dbg = struct();
dbg.fs      = fs;
dbg.modeSel = "default_fallback";
out.debug = dbg;
end


