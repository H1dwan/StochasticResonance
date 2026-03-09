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
T  = 2000;             % 信号长度 s (缩短时长以加快演示速度，实际可设为1000)
N  = fs * T;
t  = (0:N-1)'/fs;

% 目标信号
f0 = 0.01;            % 目标频率
A0 = 0.10;            % 信号幅值
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
max_test_layers = 7;  % 测试的最大层数 (例如测试 1~5 层)
num_repeats     = 1; % 每个层数重复运行次数
search_agents   = 20; % PSO 种群规模 (适当减小以节省时间)
max_iter        = 50; % PSO 迭代次数

% 参数范围 (CBSR: a, b)
lb = [0.01, 0.01];
ub = [10.00, 10.00];
dim = 2;

% 数据记录容器
record_time_avg = zeros(max_test_layers, 1);  % 每个层数的平均总耗时
record_snr_avg  = zeros(max_test_layers, 1);  % 每个层数的平均输出SNR

%% 3. 逐层优化与计时循环
fprintf('\n>>> 开始逐层优化与计时测试（每个层数重复 %d 次）...\n', num_repeats);
fprintf('| Layer | Avg Time (s) | Avg Output SNR (dB) |\n');
fprintf('|-------|--------------|---------------------|\n');

for layer = 1:max_test_layers
    snr_runs = zeros(num_repeats, 1);
    time_runs = zeros(num_repeats, 1);
    
    fprintf('  Layer %d: ', layer);
    
    for rep = 1:num_repeats
        % 每次重复都从同一初始输入开始，运行到指定层数
        current_input = clean_sig;
        current_noise = noise_seq;
        
        t_start = tic;
        best_score = NaN;
        
        for k = 1:layer
            fobj = @(params) CostFunction_CBSR(params, current_input, current_noise, fs, f0);
            [best_score, best_pos, ~] = PSO(search_agents, max_iter, lb, ub, dim, fobj);
            
            a_opt = best_pos(1);
            b_opt = best_pos(2);
            drift_func = @(x) CBSR_Dynamics(x, a_opt, b_opt);
            x_out = RK4Solver2(drift_func, current_input, current_noise, fs);
            
            current_input = x_out;
            current_noise = zeros(N, 1); % 后续层无新噪声
        end
        
        time_runs(rep) = toc(t_start);
        snr_runs(rep) = -best_score;
        
        % 进度显示：每10次打印一次，并在最后一次换行
        if mod(rep, 10) == 0 || rep == num_repeats
            fprintf('%d/%d ', rep, num_repeats);
            if rep == num_repeats
                fprintf('\n');
            end
        end
    end
    
    record_time_avg(layer) = mean(time_runs);
    record_snr_avg(layer) = mean(snr_runs);
    
    fprintf('|   %d   |   %8.2f   |       %8.4f      |\n', ...
        layer, record_time_avg(layer), record_snr_avg(layer));
end

fprintf('|-------|--------------|---------------------|\n');

%% 4. 结果可视化分析
figure('Color', 'w', 'Position', [100, 100, 900, 500]);

% 利用 yyaxis 绘制双轴图
yyaxis left
p1 = plot(1:max_test_layers, record_snr_avg, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('输出信噪比 SNR (dB)', 'FontSize', 12);
ylim([min(record_snr_avg)-2, max(record_snr_avg)+2]);
xlabel('级联层数 (Number of Layers)', 'FontSize', 12);

yyaxis right
p2 = plot(1:max_test_layers, record_time_avg, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('平均运行时间 (Seconds)', 'FontSize', 12);
ylim([0, max(record_time_avg)*1.2]);

% 美化图表
title(sprintf('级联SR系统：性能提升 vs 计算成本（每层重复 %d 次取平均）', num_repeats), 'FontSize', 14);
grid on;
xticks(1:max_test_layers);
legend([p1, p2], {'Average Output SNR', 'Average Runtime'}, 'Location', 'northwest');

% 在图上标注具体数值
for i = 1:max_test_layers
    yyaxis left
    text(i, record_snr_avg(i)+0.2, sprintf('%.1fdB', record_snr_avg(i)), ...
        'HorizontalAlignment', 'center', 'Color', 'b', 'FontSize', 10);
    
    yyaxis right
    text(i, record_time_avg(i)-max(record_time_avg)*0.05, sprintf('%.1fs', record_time_avg(i)), ...
        'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 10);
end

%% 5. 保存结果数据
results.input.A0 = A0;
results.input.f0 = f0;
results.input.fs = fs;
results.input.T  = T;
results.input.D = D;
results.input.num_repeats = num_repeats;
results.input.clean_sig = clean_sig;
results.input.noise_seq = noise_seq;
results.output.record_snr_avg = record_snr_avg;
results.output.record_time_avg = record_time_avg;

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