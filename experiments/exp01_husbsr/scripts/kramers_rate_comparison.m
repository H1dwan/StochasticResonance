% % =========================================================================
% % Description: 基于 Kramers 理论的逃逸率对比实验
% %              比较 HSUBSR、UBSR 和 PLBSR:
% %              1) 势函数得到的理论 Kramers 逃逸率 r_K(D)
% %              2) SDE 仿真得到的实际逃逸率 r_sim(D)
% %
% % Author: LiuShuang
% % Created: 2025-12-05
% % Last Modified: 2026-01-14
% %
% % Usage:
% %   修改 rate_sim_hsubsr 为 0.001185000000000
% %   直接运行本脚本即可。
% % =========================================================================
clc; clear; close all;

% 引入绘图风格设置 (如果可用)
if exist('SetThesisDefaultStyle', 'file')
    SetThesisDefaultStyle();
end

%% 1. 全局仿真参数设置 =================================================
% ---------------------------------------------------------------------
% 物理参数
% ---------------------------------------------------------------------
xm_target = 1.0;    % 目标势阱位置
dU_target = 0.25;   % 目标势垒高度
A = 0.1;            % 输入信号振幅
f = 0.01;           % 输入信号频率 (Hz)

% ---------------------------------------------------------------------
% 仿真参数
% ---------------------------------------------------------------------
fs              = 200;                % 采样率 (Hz)
t_total         = 1000;               % 单次仿真总时长 (s)，需保证足够发生跃迁
transient_ratio = 0.1;                % 瞬态截断比例
n_realizations  = 2;                  % 每个噪声强度下的平均次数 (增加以提高精度)

% 噪声强度扫描范围
D_list = 0.0:0.01:1.50;
num_D  = length(D_list);

%% 2. 模型参数校准 (确保所有模型具有相同的 xm 和 dU) ===================
fprintf('正在校准模型参数...\n');

% (1) HSUBSR 参数校准
shape_factor = 100; % 形状因子，控制势阱宽窄
[a_hs, b_hs, k1_hs, k2_hs] = CalibrateHSUBSR(xm_target, dU_target, shape_factor);

% (2) UBSR 参数校准
% UBSR: xm = sqrt(a/b), dU = a^2/(4b)
% 解得: a = 4*dU/xm, b = 4*dU/xm^2 (仅适用标准UBSR, 此处根据您的UBSR_Potential调整)
% 假设 UBSR_Potential 使用标准参数定义
a_ubsr = 4 * dU_target / (xm_target^2); 
b_ubsr = 4 * dU_target / (xm_target^4); 

% (3) PLBSR 参数校准
% PLBSR: 势阱位置 L0, 势垒高度 U0 (直接对应)
u_pl = dU_target;
l_pl = xm_target;

%% 4. 主循环：理论计算与蒙特卡洛仿真 ====================================
rate_sim_hs = zeros(num_D, 1);
rate_sim_ub = zeros(num_D, 1);
rate_sim_pl = zeros(num_D, 1);

% 定义动力学方程句柄
drift_hs = @(x) HSUBSR_Dynamics(x, a_hs, b_hs, k1_hs, k2_hs);
drift_ub = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
drift_pl = @(x) PLBSR_Dynamics(x, u_pl, l_pl);

% 判别阈值 (用于统计跃迁)
x_thr = 0.2 * xm_target; 

fprintf('\n开始仿真 (Total D points: %d)...\n', num_D);

[~, theo_snr_hs, rate_theo_hs] = HSUBSR_SNR(xm_target, dU_target, shape_factor, f, A, D_list);
[~, theo_snr_ub, rate_theo_ub] = UBSR_SNR(a_ubsr, b_ubsr, f, A, D_list);
[~, theo_snr_pl, rate_theo_pl] = PLBSR_SNR(u_pl, l_pl, f, A, D_list);

figure;
plot(D_list, theo_snr_hs, 'r-o', D_list, theo_snr_ub, 'b-s', D_list, theo_snr_pl, 'g-d', 'LineWidth', 1.5);
legend({'HSUBSR', 'UBSR', 'PLBSR'}, 'Location', 'Best');
xlabel('Noise Intensity D');
ylabel('Theoretical SNR');
title('Theoretical SNR vs Noise Intensity');

for i = 1:num_D
    D = D_list(i);
    % --- B. 仿真值统计 (蒙特卡洛) ---
    % 使用 parfor 加速 (如果并行工具箱可用，否则改为 for)
    % 这里为了演示清晰使用普通 for 循环，实际运行时建议开启并行
    % HSUBSR
    rate_sim_hs(i) = MonteCarloEscapeRate(drift_hs, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
    % UBSR
    rate_sim_ub(i) = MonteCarloEscapeRate(drift_ub, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
    % PLBSR
    rate_sim_pl(i) = MonteCarloEscapeRate(drift_pl, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
    fprintf('D = %.3f | HS_Err: %.1f%% | UB_Err: %.1f%% | PL_Err: %.1f%%\n', ...
        D, ...
        abs(rate_sim_hs(i)-rate_theo_hs(i))/rate_theo_hs(i)*100, ...
        abs(rate_sim_ub(i)-rate_theo_ub(i))/rate_theo_ub(i)*100, ...
        abs(rate_sim_pl(i)-rate_theo_pl(i))/rate_theo_pl(i)*100);
end

%% 5. 结果可视化 =======================================================
if exist('CreateThesisFigure', 'file')
    fig = CreateThesisFigure(16, 8);
else
    fig = figure('Position', [100, 100, 1000, 500], 'Color', 'w');
end

% 计算相对误差
err_hs = abs(rate_sim_hs - rate_theo_hs) ./ rate_theo_hs;
err_ub = abs(rate_sim_ub - rate_theo_ub) ./ rate_theo_ub;
err_pl = abs(rate_sim_pl - rate_theo_pl) ./ rate_theo_pl;

% 子图1：Kramers Rate (Log Scale)
subplot(1, 2, 1);
semilogy(D_list, rate_theo_hs, 'k-', 'LineWidth', 1.5); hold on;
semilogy(D_list, rate_sim_hs,  'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
semilogy(D_list, rate_theo_ub, 'b--', 'LineWidth', 1.2); 
semilogy(D_list, rate_sim_ub,  'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
semilogy(D_list, rate_theo_pl, 'g-.', 'LineWidth', 1.2);
semilogy(D_list, rate_sim_pl,  'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
grid on;
xlabel('Noise Intensity D');
ylabel('Kramers Escape Rate (log scale)');
title('Theoretical vs Simulated Rate');
legend({'HS-Theory', 'HS-Sim', 'UB-Theory', 'UB-Sim', 'PL-Theory', 'PL-Sim'}, ...
    'Location', 'SouthEast', 'FontSize', 9);

% 子图2：相对误差对比
subplot(1, 2, 2);
plot(D_list, err_hs, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
plot(D_list, err_ub, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
plot(D_list, err_pl, 'g-d', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
grid on;
xlabel('Noise Intensity D');
ylabel('Relative Error (|r_{sim} - r_{th}| / r_{th})');
title('Deviation from Theory');
legend({'HSUBSR', 'UBSR', 'PLBSR'}, 'Location', 'Best');
ylim([0, 1]); % 限制误差显示范围以便观察

%% ====================================================================
%  Local Helper Functions
% ====================================================================
function rate_mean = MonteCarloEscapeRate(drift_func, fs, t_total, d_noise, ...
    x_thr, n_realizations, transient_ratio)
% MonteCarloEscapeRate  基于多次仿真估计逃逸率
%
% 输入:
%   drift_func       - 漂移函数句柄 f(x)
%   fs               - 采样率 (Hz)
%   t_total          - 单次仿真总时长 (s)
%   d_noise          - 噪声强度 D
%   x_thr            - 判别左右势阱的阈值 |x| > x_thr 视为在某个势阱内
%   n_realizations   - 重复仿真次数
%   transient_ratio  - 瞬态比例 (0~1)，前期数据丢弃不用统计
%
% 输出:
%   rate_mean        - 该模型在该 D 下的平均逃逸率估计

h_step   = 1 / fs;
n_steps  = round(t_total * fs);
n_points = n_steps + 1;

clean_signal = zeros(n_points, 1);     % 无输入信号 s(t) = 0
idx0         = round(transient_ratio * n_points) + 1;

rate_list = zeros(n_realizations, 1);

for i_real = 1:n_realizations
    % 生成高斯白噪声序列 noise_seq，使得:
    % dx_stoch ≈ noise_seq(i) * h_step = sqrt(2D h) * N(0,1)
    noise_seq = sqrt(2 * d_noise / h_step) * randn(n_points, 1);
    
    x_traj = RK4Solver2(drift_func, clean_signal, noise_seq, fs);
    
    % 丢弃瞬态，仅用稳态部分估计逃逸率
    x_stead = x_traj(idx0:end);
    
    rate_list(i_real) = EstimateEscapeRate(x_stead, fs, x_thr);
end

rate_mean = mean(rate_list);
end

function rate_est = EstimateEscapeRate(x_traj, fs, x_thr)
% EstimateEscapeRate  根据势阱跃迁次数估计逃逸率
%
% 算法:
%   1. 定义三类区域:
%      左阱: x <= -x_thr
%      右阱: x >=  x_thr
%      中间: 其他 (靠近势垒，不计入势阱驻留)
%   2. 沿时间扫描区域标签，统计从左阱到右阱或从右阱到左阱的跃迁次数。
%   3. 逃逸率估计为:
%      r_hat = (跃迁总次数) / T_obs
%
%   对称双稳系统的两阱跃迁率相同，且该定义下的“单位时间跃迁总数”
%   即为单方向的 Kramers 逃逸率 r。

n_points = numel(x_traj);
region   = zeros(n_points, 1);

% 区域标签: -1 = 左阱, +1 = 右阱, 0 = 中间区
region(x_traj <= -x_thr) = -1;
region(x_traj >=  x_thr) =  1;

last_region   = 0;
trans_count   = 0;

for i = 1:n_points
    current_region = region(i);
    
    if current_region == 0
        % 在中间过渡区，不更新跃迁计数，只是停留
        continue;
    end
    
    if last_region == 0
        % 第一次进入某个势阱
        last_region = current_region;
    elseif current_region ~= last_region
        % 从左阱跳到右阱，或从右阱跳到左阱，计一次跃迁
        trans_count = trans_count + 1;
        last_region = current_region;
    end
end

t_obs   = (n_points - 1) / fs;
rate_est = trans_count / t_obs;
end