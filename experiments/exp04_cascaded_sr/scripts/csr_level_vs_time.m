% =========================================================================
% Description: 级联SR系统性能与计算耗时分析 (Performance vs. Time Cost)
%
% Author: LiuShuang
% Created: 2026-01-09
% Last Modified: 2026-01-09
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 仿真信号与环境设置
fs = 5;              % 采样频率 Hz
T  = 1000;             % 信号长度 s (缩短时长以加快演示速度，实际可设为1000)
N  = fs * T;
t  = (0:N-1)'/fs;

% 目标信号
f0 = 0.01;            % 目标频率
A0 = 0.05;            % 信号幅值
clean_sig = A0 * sin(2*pi*f0*t);

% 背景噪声
D  = 1.0;             % 噪声强度
h  = 1/fs;
noise_seq = sqrt(2*D/h) * randn(N, 1);

input_sig = clean_sig + noise_seq;
snr_in = SNRo2(input_sig, fs, f0);

fprintf('>>> 仿真环境配置:\n');
fprintf('    采样率: %d Hz, 时长: %d s\n', fs, T);
fprintf('    输入 SNR: %.4f dB\n', snr_in);

%% 2. 优化与计时参数
max_test_layers = 6;  % 测试的最大层数 (例如测试 1~5 层)
search_agents   = 10; % PSO 种群规模 (适当减小以节省时间)
max_iter        = 20; % PSO 迭代次数

% 参数范围 (CBSR: a, b)
lb = [0.01, 0.01];
ub = [10.00, 10.00];
dim = 2;

% 数据记录容器
record_time_step = zeros(max_test_layers, 1); % 每一层的单步耗时
record_time_cum  = zeros(max_test_layers, 1); % 累计耗时 (N层系统总耗时)
record_snr       = zeros(max_test_layers, 1); % 每一层的输出SNR

% 初始化逐级循环变量
current_input = clean_sig;
current_noise = noise_seq;
total_elapsed_time = 0;

%% 3. 逐层优化与计时循环
fprintf('\n>>> 开始逐层优化与计时测试...\n');
fprintf('| Layer | Step Time (s) | Total Time (s) | Output SNR (dB) |\n');
fprintf('|-------|---------------|----------------|-----------------|\n');

for layer = 1:max_test_layers
    % --- 计时开始 ---
    t_start = tic;
    
    % 1. 定义适应度函数 (闭包捕获当前输入)
    fobj = @(params) CostFunction_CBSR(params, current_input, current_noise, fs, f0);
    
    % 2. 执行 PSO 优化
    % 注意：PSO 内部有随机性，但计算量相对固定
    [best_score, best_pos, ~] = PSO(search_agents, max_iter, lb, ub, dim, fobj);
    
    % 3. 使用最优参数生成输出 (用于下一层)
    a_opt = best_pos(1);
    b_opt = best_pos(2);
    drift_func = @(x) CBSR_Dynamics(x, a_opt, b_opt);
    x_out = RK4Solver2(drift_func, current_input, current_noise, fs);
    
    % --- 计时结束 ---
    t_step = toc(t_start);
    
    % --- 更新数据 ---
    total_elapsed_time = total_elapsed_time + t_step;
    
    record_time_step(layer) = t_step;
    record_time_cum(layer)  = total_elapsed_time;
    record_snr(layer)       = -best_score; % 还原 SNR
    
    % 打印当前层结果
    fprintf('|   %d   |    %7.2f    |    %7.2f     |     %7.4f     |\n', ...
        layer, t_step, total_elapsed_time, record_snr(layer));
    
    % 更新下一层输入 (级联操作)
    current_input = x_out;
    current_noise = zeros(N, 1); % 后续层无新噪声
end

fprintf('|-------|---------------|----------------|-----------------|\n');

%% 4. 结果可视化分析
figure('Color', 'w', 'Position', [100, 100, 900, 500]);

% 利用 yyaxis 绘制双轴图
yyaxis left
p1 = plot(1:max_test_layers, record_snr, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('输出信噪比 SNR (dB)', 'FontSize', 12);
ylim([min(record_snr)-2, max(record_snr)+2]);
xlabel('级联层数 (Number of Layers)', 'FontSize', 12);

yyaxis right
p2 = plot(1:max_test_layers, record_time_cum, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('累计计算耗时 (Seconds)', 'FontSize', 12);
ylim([0, max(record_time_cum)*1.2]);

% 美化图表
title('级联SR系统：性能提升 vs 计算成本', 'FontSize', 14);
grid on;
xticks(1:max_test_layers);
legend([p1, p2], {'Output SNR', 'Total Computation Time'}, 'Location', 'northwest');

% 在图上标注具体数值
for i = 1:max_test_layers
    yyaxis left
    text(i, record_snr(i)+0.2, sprintf('%.1fdB', record_snr(i)), ...
        'HorizontalAlignment', 'center', 'Color', 'b', 'FontSize', 10);
    
    yyaxis right
    text(i, record_time_cum(i)-max(record_time_cum)*0.05, sprintf('%.1fs', record_time_cum(i)), ...
        'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 10);
end

%% 5. 保存结果数据
results.input.A0 = A0;
results.input.f0 = f0;
results.input.fs = fs;
results.input.T  = T;
results.input.D = D;
results.input.clean_sig = clean_sig;
results.input.noise_seq = noise_seq;
results.output.record_snr = record_snr;
results.output.record_time_cum = record_time_cum;

%% 辅助函数 (Cost Function)
function fitness = CostFunction_CBSR(params, sig_in, noise_in, fs, f0)
a = params(1);
b = params(2);
drift = @(x) CBSR_Dynamics(x, a, b);
try
    x_out = RK4Solver2(drift, sig_in, noise_in, fs);
    N = length(x_out);
    x_steady = x_out(round(0.1*N):end); % 去除瞬态
    val = SNRo2(x_steady, fs, f0);
    fitness = -val;
    if isnan(fitness), fitness = 1e5; end
catch
    fitness = 1e5;
end
end