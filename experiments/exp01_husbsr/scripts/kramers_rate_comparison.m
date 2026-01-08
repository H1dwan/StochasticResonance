% =========================================================================
% Description: 基于 Kramers 理论的逃逸率对比实验
%              比较 HSUBSR、UBSR 和 PLBSR:
%              1) 势函数得到的理论 Kramers 逃逸率 r_K(D)
%              2) SDE 仿真得到的实际逃逸率 r_sim(D)
%
% Author: LiuShuang
% Created: 2025-12-05
% Last Modified: 2025-12-22
%
% Usage:
%   修改 rate_sim_hsubsr 为 0.001185000000000
%   直接运行本脚本即可。
% =========================================================================

clc; clear; close all;
rng(42);  % 固定随机种子，保证可重复性

%% 1. 全局仿真参数与模型参数 ==========================================

% -------- 常量设置 ------------------------------------------------------
T_TOTAL          = 2000;      % 总仿真时间 (s)，保证足够多的跃迁
FS               = 200;       % 采样率 (Hz)
TRANSIENT_RATIO  = 0.2;       % 前 20% 作为瞬态丢弃
N_REALIZATIONS   = 2;         % 每个 D 的重复次数，降低统计方差

D_LIST           = 0.05:0.05:0.50;  % 噪声强度列表
NUM_D            = numel(D_LIST);

% -------- 固定的势函数结构参数 ----------------------------------------
xm = 1;
dU = 0.25;

% -------- HSUBSR 势函数参数 ------------------------------------------
shape_factor = 6;
[a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr] = CalibrateHSUBSR(xm, dU, shape_factor);

% -------- UBSR 势函数参数 --------------------------------------------
a_ubsr     = 1.0;
b_ubsr     = 1.0;

% -------- PLBSR 势函数参数 -------------------------------------------
u_plbsr    = 0.25;
l_plbsr    = 1;

%% 2. 计算模型的 Kramers 理论参数 (ΔU, U'' 等) =====================

[xm_hsubsr, xb_hsubsr, delta_u_hsubsr, uxx_m_hsubsr, uxx_b_hsubsr] = ...
    ComputeKramersParamsHSUBSR(a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);

[xm_ubsr, xb_ubsr, delta_u_ubsr, uxx_m_ubsr, uxx_b_ubsr] = ...
    ComputeKramersParamsUBSR(a_ubsr, b_ubsr);

% [xm_pl, xb_pl, delta_u_pl, uxx_m_pl, uxx_b_pl] = ...
%     ComputeKramersParamsPLBSR(u_plbsr, l_plbsr);

% Kramers 理论中的前因子 (与 D 无关)
prefactor_hsubsr = sqrt(abs(uxx_b_hsubsr) * uxx_m_hsubsr) / (2 * pi);
% prefactor_ubsr = sqrt(abs(uxx_b_ubsr) * uxx_m_ubsr) / (2 * pi);
prefactor_ubsr = a_ubsr^2 / (4*b_ubsr*sqrt(a_ubsr/b_ubsr) * (sqrt(2*a_ubsr/b_ubsr) - sqrt(a_ubsr/b_ubsr)));
% prefactor_plbsr = sqrt(abs(uxx_b_pl) * uxx_m_pl) / (2 * pi);
prefactor_plbsr = u_plbsr^2 / (4*l_plbsr^2);

% 定义阈值，用于区分“左阱 / 右阱”区域，避免在势垒附近误计多次跃迁
x_thr_hsubsr = 0.5 * xm_hsubsr;
x_thr_ubsr = 0.5 * xm;
x_thr_plbsr  = 0.5 * xm;

%% 3. 不同 D 下的理论 Kramers 逃逸率曲线 ==============================

rate_theory_hsubsr = zeros(NUM_D, 1);
rate_theory_ubsr = zeros(NUM_D, 1);
rate_theory_plbsr  = zeros(NUM_D, 1);

for i_d = 1:NUM_D
    d_val = D_LIST(i_d);
    rate_theory_hsubsr(i_d) = prefactor_hsubsr * exp(-delta_u_hsubsr / d_val);
    rate_theory_ubsr(i_d) = prefactor_ubsr * exp(-dU / d_val);
    rate_theory_plbsr(i_d) = prefactor_plbsr * (1 / (d_val * sinh(u_plbsr/2/d_val)^2));
end

%% 4. 基于 SDE 仿真的逃逸率估计 (不用 ZCR，显式数跃迁次数) =============

rate_sim_hsubsr = zeros(NUM_D, 1);
rate_sim_ubsr = zeros(NUM_D, 1);
rate_sim_plbsr  = zeros(NUM_D, 1);

drift_hsubsr_func = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);
drift_ubsr_func = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
drift_plbsr_func = @(x) PLBSR_Dynamics(x, u_plbsr, l_plbsr);
for i_d = 1:NUM_D
    d_val = D_LIST(i_d);
    
    rate_sim_hsubsr(i_d) = MonteCarloEscapeRate( ...
        drift_hsubsr_func, FS, T_TOTAL, d_val, x_thr_hsubsr, ...
        N_REALIZATIONS, TRANSIENT_RATIO);
    
    rate_sim_ubsr(i_d) = MonteCarloEscapeRate( ...
        drift_ubsr_func, FS, T_TOTAL, d_val, x_thr_ubsr, ...
        N_REALIZATIONS, TRANSIENT_RATIO);
    
    rate_sim_plbsr(i_d) = MonteCarloEscapeRate( ...
        drift_plbsr_func, FS, T_TOTAL, d_val, x_thr_plbsr, ...
        N_REALIZATIONS, TRANSIENT_RATIO);
    
    fprintf('D = %.3f:  HSUBSR r_sim = %.4e,  UBSR r_sim = %.4e,  PLBSR r_sim = %.4e\n', ...
        d_val, rate_sim_hsubsr(i_d), rate_sim_ubsr(i_d), rate_sim_plbsr(i_d));
end

%% 5. 计算相对误差 |r_sim - r_K| / r_K ==================================

rel_err_hsubsr = abs(rate_sim_hsubsr - rate_theory_hsubsr) ./ rate_theory_hsubsr;
rel_err_ubsr = abs(rate_sim_ubsr - rate_theory_ubsr) ./ rate_theory_ubsr;
rel_err_plbsr = abs(rate_sim_plbsr - rate_theory_plbsr) ./ rate_theory_plbsr;

%% 6. 绘图：Kramers 理论 vs 仿真 & 相对误差 ============================
SetThesisDefaultStyle();
% fig = tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
fig = CreateThesisFigure(16,6);

% ---- 左图：逃逸率-噪声强度 (半对数) -------------------------------
subplot(1,2,1);
semilogy(D_LIST, rate_theory_hsubsr, '-', 'LineWidth', 1.5, 'Color', [0.9020, 0.2941, 0.2078]); hold on;
semilogy(D_LIST, rate_sim_hsubsr,   'o', 'LineWidth', 1.2, 'Color', [0.9020, 0.2941, 0.2078]);
semilogy(D_LIST, rate_theory_ubsr, '-', 'LineWidth', 1.5, 'Color', [0.2510, 0.3686, 0.6353]);
semilogy(D_LIST, rate_sim_ubsr,   's', 'LineWidth', 1.2, 'Color', [0.2510, 0.3686, 0.6353]);
semilogy(D_LIST, rate_theory_plbsr,  '-', 'LineWidth', 1.5, 'Color', [0.0000, 0.6196, 0.4510]);
semilogy(D_LIST, rate_sim_plbsr,    'd',  'LineWidth', 1.2, 'Color', [0.0000, 0.6196, 0.4510]);
xlabel('Noise Intensity D');
ylabel('Kramers Rate');
title('Comparison of Kramers Rate');
legend({'HSUBSR theory','HSUBSR sim', ...
    'UBSR theory','UBSR sim', ...
    'PLBSR  theory','PLBSR  sim'}, 'Location', 'best');
grid on;

% ---- 右图：相对误差对比 ---------------------------------------------
subplot(1,2,2);
plot(D_LIST, rel_err_hsubsr, 'o-', 'LineWidth', 1.5); hold on;
plot(D_LIST, rel_err_ubsr, 's--', 'LineWidth', 1.5);
plot(D_LIST, rel_err_plbsr,  'd-.', 'LineWidth', 1.5);
xlabel('Noise Intensity D');
ylabel('Relative Error');
title('Fitting Error');
legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'best');
grid on;


%% 辅助函数部分
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

function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
    ComputeKramersParamsHSUBSR(a, b, k1, k2)
% ComputeKramersParamsHSUBSR  计算 HSUBSR 势函数的 Kramers 参数
%
% 输出:
%   x_m     - 右侧势阱最小值位置
%   x_b     - 中心势垒位置 (对称模型中为 0)
%   delta_u - 势垒高度 ΔU = U(x_b) - U(x_m)
%   uxx_m   - U''(x_m)
%   uxx_b   - U''(x_b) (一般为负值)
U_fun = @(x) HSUBSR_Potential(x, a, b, k1, k2);

% 对称双稳势，势垒在 0，右侧势阱在 x > 0
x_search_right = max(3 * 1, 3);   % 搜索区间上限，可根据参数调整
x_m = fminbnd(U_fun, 0.1, x_search_right);
x_b = 0.0;

% 势垒高度
delta_u = U_fun(x_b) - U_fun(x_m);

% 数值二阶导 (有限差分)
h_fd  = 1e-3;
uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
end

function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
    ComputeKramersParamsUBSR(a, b)
% ComputeKramersParamsUBSR  计算 UBSR 势函数的 Kramers 参数
%
% 注意:
%   需要 MATLAB 路径中存在 UBSR_Potential.m
U_fun    = @(x) UBSR_Potential(x, 'a', a, 'b', b);
threshold = sqrt(a / b);

% 右侧势阱在 x > 0，约在 threshold 附近
x_m = fminbnd(U_fun, 0, 2 * threshold);
x_b = 0.0;   % 对称模型中心势垒

delta_u = U_fun(x_b) - U_fun(x_m);

h_fd  = 1e-3;
uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
end

% function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
%     ComputeKramersParamsPLBSR(a, b, c)
% % ComputeKramersParamsPLBSR  计算 PLBSR 势函数的“有效” Kramers 参数
% %
% % 说明:
% %   PLBSR 为分段线性势函数，在 x = 0, ±b2 处二阶导数不存在，
% %   不严格满足 Kramers 理论的光滑假设。这里采用有限差分在
% %   极值点附近给出一个“有效曲率”，用作对比。

% U_fun    = @(x) PLBSR_Potential(x, 'a2', a, 'b2', b, 'c2', c);
% % 右侧势阱理论上在 x = b
% x_m = b;
% x_b = 0.0;          % 中心势垒

% delta_u = U_fun(x_b) - U_fun(x_m);
% % 由于为分段线性，二阶导数的有限差分近似依赖于步长 h_fd，
% % 这里只取一个相对“粗”的 h_fd 作为有效曲率估计。
% h_fd  = 1e-2;
% uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
% uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
% end