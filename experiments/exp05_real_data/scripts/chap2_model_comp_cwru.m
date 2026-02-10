% =========================================================================
% Description: 对比UBSR PLBSR HSUBSR三种模型在CWRU轴承数据上的性能
%
% Author: LiuShuang
% Created: 2026-02-10
% Last Modified: 2026-02-10
%
% Usage: 保证路径下存在CWRU数据文件
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;
seed = 5;

% 1. 加载数据 ==============================================================
data_filename = '105.mat';
S = load(data_filename);  % 加载到结构体

% 查找变量名
var_names = fieldnames(S);
de_idx = find(contains(var_names, '_DE_time'));
rpm_idx = find(contains(var_names, 'RPM'));

if isempty(de_idx) || isempty(rpm_idx)
    error('找不到所需的变量');
end

% 提取数据
raw_sig = S.(var_names{de_idx(1)});
rpm = S.(var_names{rpm_idx(1)});

% 计算轴承故障频率
[freqs, desc] = BearingHz(rpm);

% 2. 信号预处理 ============================================================
% 截取数据片段进行分析（避免计算量过大）
fs_raw = 12000;      % CWRU标准采样率 12kHz
N_sample = 8192; 
N = length(raw_sig);
t = (0:N-1)/fs_raw;
sig_segment = raw_sig(1:N_sample);

% 包络提取，SR作为低通滤波器，直接处理高频载波效果差，需先提取故障包络
sig_env = abs(hilbert(sig_segment));

% 2.2 去直流
sig_ac = sig_env - mean(sig_env);

% 2.3 幅值尺度变换 (Amplitude Rescaling)
% z-score 标准化
sig_mean = mean(sig_ac);
sig_std = std(sig_ac);
sig_processed = (sig_ac - sig_mean) / sig_std;

% 2.4 频率尺度变换 (Frequency Rescaling)
% 将真实采样频率映射到固定的5Hz 
f_fault_real = freqs.BPFI;      % 目标故障频率
R = fs_raw / 5;     % 尺度变换系数 R
f_target_sr = f_fault_real / R;     % 变换后的目标频率

% 更新数值积分步长 h
% 原始步长 h0 = 1/fs. 变换后系统"感知"到的步长需放大 R 倍
% 理论依据: d x / d (t/R) = ... => 相当于时间变慢，频率变低
h_sr = (1/fs_raw) * R;

%% 3. 自适应 PSO 寻优 (Adaptive Optimization) ==============================

% PSO 参数设置
dim = 2;              % 优化变量维度 [xm, dU]
lb = [0.1, 0.1];      % 参数下界
ub = [5.0, 5.0];      % 参数上界

dim_hs = 3;           % HSUBSR 多一个 shape 参数
lb_hs = [0.1, 0.1, 1.01];
ub_hs = [5.0, 5.0, 50];

SearchAgents_no = 20; % 种群规模
Max_iter = 50;        % 迭代次数

% 定义目标函数: 输入[a,b], 输出 -SNR (因为PSO默认求极小值)
fobj_ubsr  = @(p) Fitness_UBSR_Fair(p, sig_processed, h_sr, f_target_sr);
fobj_plbsr = @(p) Fitness_PLBSR_Fair(p, sig_processed, h_sr, f_target_sr);
fobj_hs    = @(p) Fitness_HSUBSR_Fair(p, sig_processed, h_sr, f_target_sr);

fprintf('\nStarting PSO Optimization...\n');

% 并行计算策略：由于PSO具有随机性，运行多次取最优
n_pso_runs = 1;

models = struct( ...
    'name', {'UBSR', 'PLBSR', 'HSUBSR'}, ...
    'dim', {dim, dim, dim_hs}, ...
    'lb',  {lb,  lb,  lb_hs}, ...
    'ub',  {ub,  ub,  ub_hs}, ...
    'fobj', {fobj_ubsr, fobj_plbsr, fobj_hs} ...
    );

results = repmat(struct('name', '', 'best_snr', -inf, 'best_params', [], ...
    'curve', [], 'runs', []), numel(models), 1);

tic;
for m = 1:numel(models)
    run_log = cell(n_pso_runs, 1);
    for i = 1:n_pso_runs
        % 每个 worker 独立运行一次完整的 PSO
        [best_score, best_pos, curve] = PSO(SearchAgents_no, Max_iter, ...
            models(m).lb, models(m).ub, models(m).dim, models(m).fobj, false, seed);
        % 记录结果 (注意: PSO返回的是 minimized value，即 -SNR)
        current_max_snr = -best_score;
        run_log{i} = struct('snr', current_max_snr, 'params', best_pos, 'curve', curve);
    end

    % 汇总当前模型的多次结果
    best_snr = -inf;
    best_params = [];
    best_curve = [];
    for i = 1:n_pso_runs
        if run_log{i}.snr > best_snr
            best_snr = run_log{i}.snr;
            best_params = run_log{i}.params;
            best_curve = run_log{i}.curve;
        end
    end

    results(m).name = models(m).name;
    results(m).best_snr = best_snr;
    results(m).best_params = best_params;
    results(m).curve = best_curve;
    results(m).runs = run_log;
end
toc;

fprintf('Optimization Finished.\n');
for m = 1:numel(results)
    fprintf('  -> %s Best SNR: %.2f dB\n', results(m).name, results(m).best_snr);
    fprintf('     Params: %s\n', mat2str(results(m).best_params, 4));
end

% 最优结果验证与可视化
for m = 1:numel(results)
    fprintf('\nValidating Model: %s\n', results(m).name);
    best_params = results(m).best_params;
    switch results(m).name
        case 'UBSR'
            x_sr = GetOutput_UBSR_Fair(best_params, sig_processed, h_sr);
        case 'PLBSR'
            x_sr = GetOutput_PLBSR_Fair(best_params, sig_processed, h_sr);
        case 'HSUBSR'
            x_sr = GetOutput_HSUBSR_Fair(best_params, sig_processed, h_sr);
        otherwise
            error('Unknown model name: %s', results(m).name);
    end
    Plot_Time_Frequency(x_sr, 5, length(x_sr));
end

%% 辅助函数
function [freqs, description] = BearingHz(Fr_rpm, D, d, Z, alpha_deg, verbose)
% BearingHz 计算滚动轴承的故障特征频率
% 
% Useage:
%   freqs = BearingHz(1797); % 使用 CWRU 默认参数 (SKF 6205-2RS)
%   freqs = BearingHz(1797, 39.04, 7.94, 9, 0); % 自定义参数
%
% Input:
%   Fr_rpm    : 轴旋转速度 (RPM), 支持向量输入
%   D         : 节径 (Pitch Diameter), 单位 mm (默认: CWRU参数 39.04mm)
%   d         : 滚动体直径 (Ball Diameter), 单位 mm (默认: CWRU参数 7.94mm)
%   Z         : 滚动体个数 (Number of Elements) (默认: 9)
%   alpha_deg : 接触角 (Contact Angle), 单位 度 (默认: 0)
%   verbose   : 是否打印结果 (true/false, 默认 true)
%
% Output:
%   freqs     : 包含各故障频率的结构体 (Hz)
%       .BPFI : 内圈故障频率
%       .BPFO : 外圈故障频率
%       .BSF  : 滚动体故障频率
%       .FTF  : 保持架(基频)频率
%       .Fr   : 转频
%
% Reference: 
%   CWRU Bearing: SKF 6205-2RS JEM SKF
%   D = 1.537 inch ≈ 39.04 mm
%   d = 0.3126 inch ≈ 7.94 mm
%   Z = 9 balls
%   alpha = 0 deg

% =========================================================================
% 1. 参数预处理与默认值设置 (针对 CWRU 数据集)
% =========================================================================
if nargin < 6, verbose = true; end
% 如果只输入了转速，默认加载 SKF 6205-2RS 参数
if nargin < 5 || (isempty(D) && isempty(d))
    % SKF 6205-2RS 参数 (CWRU标准)
    D = 39.04; 
    d = 7.94;  
    Z = 9;     
    alpha_deg = 0; 
    if verbose && nargin < 2
        fprintf('Info: Using default CWRU (SKF 6205-2RS) parameters.\n');
    end
end

% 角度转弧度
alpha = alpha_deg * pi / 180;

% 转速 RPM -> Hz
fr = Fr_rpm / 60; 

% 几何比率 (简化公式书写)
ratio = (d / D) * cos(alpha);

% =========================================================================
% 2. 频率计算 (Hz)
% =========================================================================

% 内圈故障频率 (Inner Race)
bpfi = (Z / 2) * (1 + ratio) .* fr;

% 外圈故障频率 (Outer Race)
bpfo = (Z / 2) * (1 - ratio) .* fr;

% 保持架故障频率 (Cage / Fundamental Train)
ftf  = (1 / 2) * (1 - ratio) .* fr;

% 滚动体故障频率 (Ball Spin)
% 公式: D/(2d) * (1 - (d/D*cos(a))^2) * fr
bsf  = (D / (2 * d)) * (1 - ratio.^2) .* fr;

% =========================================================================
% 3. 结果封装与显示
% =========================================================================
freqs.Fr   = fr;
freqs.BPFI = bpfi;
freqs.BPFO = bpfo;
freqs.BSF  = bsf;
freqs.FTF  = ftf;

% 生成描述文本（可选）
description = sprintf(['Rotation Frequencies (Fr = %.2f Hz):\n' ...
                       '  Inner Race (BPFI): %.2f Hz (%.2f X)\n' ...
                       '  Outer Race (BPFO): %.2f Hz (%.2f X)\n' ...
                       '  Ball Spin  (BSF) : %.2f Hz (%.2f X)\n' ...
                       '  Cage       (FTF) : %.2f Hz (%.2f X)'], ...
                       mean(fr), ...
                       mean(bpfi), mean(bpfi)/mean(fr), ...
                       mean(bpfo), mean(bpfo)/mean(fr), ...
                       mean(bsf),  mean(bsf)/mean(fr), ...
                       mean(ftf),  mean(ftf)/mean(fr));

if verbose
    fprintf('--------------------------------------------------\n');
    fprintf('%s\n', description);
    fprintf('--------------------------------------------------\n');
end

end

% --- 1. UBSR 映射与适应度 ---
function [a, b] = MapParams_UBSR(xm, dU)
    % 物理推导: U(x) = -a/2 x^2 + b/4 x^4
    % xm^2 = a/b; dU = a^2/(4b)
    % 解得: a = 4*dU / xm^2; b = 4*dU / xm^4
    a = 4 * dU / (xm^2);
    b = 4 * dU / (xm^4);
end

function val = Fitness_UBSR_Fair(p, s, h, f0)
    [a, b] = MapParams_UBSR(p(1), p(2)); % p=[xm, dU]
    drift = @(x) UBSR_Dynamics(x, a, b);
    x = RK4Solver(drift, s, h);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_UBSR_Fair(p, s, h)
    [a, b] = MapParams_UBSR(p(1), p(2));
    drift = @(x) UBSR_Dynamics(x, a, b);
    x = RK4Solver(drift, s, h);
end

% --- 2. PLBSR 映射与适应度 ---
function [U0, L0] = MapParams_PLBSR(xm, dU)
    % PLBSR 参数直接对应: L0是阱宽(xm), U0是垒高(dU)
    L0 = xm;
    U0 = dU;
end

function val = Fitness_PLBSR_Fair(p, s, h, f0)
    [U0, L0] = MapParams_PLBSR(p(1), p(2)); % p=[xm, dU]
    drift = @(x) PLBSR_Dynamics(x, U0, L0);
    x = RK4Solver(drift, s, h);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_PLBSR_Fair(p, s, h)
    [U0, L0] = MapParams_PLBSR(p(1), p(2));
    drift = @(x) PLBSR_Dynamics(x, U0, L0);
    x = RK4Solver(drift, s, h);
end

% --- 3. HSUBSR 映射与适应度 ---
function val = Fitness_HSUBSR_Fair(p, s, h, f0)
    % p = [xm, dU, shape]
    % 调用已有的校准函数
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    x = RK4Solver(drift, s, h);
    x = x(round(0.1*end):end);
    val = -SNRo2(x, 1/h, f0);
end

function x = GetOutput_HSUBSR_Fair(p, s, h)
    [a, b, k1, k2] = CalibrateHSUBSR(p(1), p(2), p(3));
    drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
    x = RK4Solver(drift, s, h);
end