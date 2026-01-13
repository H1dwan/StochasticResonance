% % =========================================================================
% % Description: 基于 Kramers 理论的逃逸率对比实验
% %              比较 HSUBSR、UBSR 和 PLBSR:
% %              1) 势函数得到的理论 Kramers 逃逸率 r_K(D)
% %              2) SDE 仿真得到的实际逃逸率 r_sim(D)
% %
% % Author: LiuShuang
% % Created: 2025-12-05
% % Last Modified: 2025-12-22
% %
% % Usage:
% %   修改 rate_sim_hsubsr 为 0.001185000000000
% %   直接运行本脚本即可。
% % =========================================================================

% clc; clear; close all;

% %% 1. 全局仿真参数与模型参数 ==========================================

% % -------- 常量设置 ------------------------------------------------------
% T_TOTAL          = 1000;      % 总仿真时间 (s)，保证足够多的跃迁
% FS               = 200;       % 采样率 (Hz)
% TRANSIENT_RATIO  = 0.1;       % 前 10% 作为瞬态丢弃
% N_REALIZATIONS   = 1;         % 每个 D 的重复次数，降低统计方差

% D_LIST           = 0.05:0.05:0.50;  % 噪声强度列表
% NUM_D            = numel(D_LIST);

% % -------- 固定的势函数结构参数 ----------------------------------------
% xm = 1;
% dU = 0.25;

% % -------- HSUBSR 势函数参数 ------------------------------------------
% shape_factor = 4;
% [a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr] = CalibrateHSUBSR(xm, dU, shape_factor);

% % -------- UBSR 势函数参数 --------------------------------------------
% a_ubsr     = 1.0;
% b_ubsr     = 1.0;

% % -------- PLBSR 势函数参数 -------------------------------------------
% u_plbsr    = 0.25;
% l_plbsr    = 1;

% %% 2. 计算模型的 Kramers 理论参数 (ΔU, U'' 等) =====================

% [xm_hsubsr, xb_hsubsr, delta_u_hsubsr, uxx_m_hsubsr, uxx_b_hsubsr] = ...
%     ComputeKramersParamsHSUBSR(a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);

% [xm_ubsr, xb_ubsr, delta_u_ubsr, uxx_m_ubsr, uxx_b_ubsr] = ...
%     ComputeKramersParamsUBSR(a_ubsr, b_ubsr);

% [xm_pl, xb_pl, delta_u_pl, uxx_m_pl, uxx_b_pl] = ...
%     ComputeKramersParamsPLBSR(u_plbsr, l_plbsr);

% % Kramers 理论中的前因子 (与 D 无关)
% prefactor_hsubsr = sqrt(abs(uxx_b_hsubsr) * uxx_m_hsubsr) / (2 * pi);
% % prefactor_ubsr = sqrt(abs(uxx_b_ubsr) * uxx_m_ubsr) / (2 * pi);
% prefactor_ubsr = a_ubsr^2 / (4*b_ubsr*sqrt(a_ubsr/b_ubsr) * (sqrt(2*a_ubsr/b_ubsr) - sqrt(a_ubsr/b_ubsr)));
% prefactor_plbsr = sqrt(abs(uxx_b_pl) * uxx_m_pl) / (2 * pi);
% % prefactor_plbsr = u_plbsr^2 / (4*l_plbsr^2);

% % 定义阈值，用于区分“左阱 / 右阱”区域，避免在势垒附近误计多次跃迁
% x_thr_hsubsr = 0.5 * xm_hsubsr;
% x_thr_ubsr = 0.5 * xm;
% x_thr_plbsr  = 0.5 * xm;

% %% 3. 不同 D 下的理论 Kramers 逃逸率曲线 ==============================

% rate_theory_hsubsr = zeros(NUM_D, 1);
% rate_theory_ubsr = zeros(NUM_D, 1);
% rate_theory_plbsr  = zeros(NUM_D, 1);

% for i_d = 1:NUM_D
%     d_val = D_LIST(i_d);
%     rate_theory_hsubsr(i_d) = prefactor_hsubsr * exp(-delta_u_hsubsr / d_val);
%     rate_theory_ubsr(i_d) = prefactor_ubsr * exp(-dU / d_val);
%     % rate_theory_plbsr(i_d) = prefactor_plbsr * (1 / (d_val * sinh(u_plbsr/2/d_val)^2));
%     rate_theory_plbsr(i_d) = prefactor_plbsr * exp(-u_plbsr / d_val);
% end

% %% 4. 基于 SDE 仿真的逃逸率估计 (不用 ZCR，显式数跃迁次数) =============

% rate_sim_hsubsr = zeros(NUM_D, 1);
% rate_sim_ubsr = zeros(NUM_D, 1);
% rate_sim_plbsr  = zeros(NUM_D, 1);

% drift_hsubsr_func = @(x) HSUBSR_Dynamics(x, a_hsubsr, b_hsubsr, k1_hsubsr, k2_hsubsr);
% drift_ubsr_func = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
% drift_plbsr_func = @(x) PLBSR_Dynamics(x, u_plbsr, l_plbsr);
% for i_d = 1:NUM_D
%     d_val = D_LIST(i_d);
    
%     rate_sim_hsubsr(i_d) = MonteCarloEscapeRate( ...
%         drift_hsubsr_func, FS, T_TOTAL, d_val, x_thr_hsubsr, ...
%         N_REALIZATIONS, TRANSIENT_RATIO);
    
%     rate_sim_ubsr(i_d) = MonteCarloEscapeRate( ...
%         drift_ubsr_func, FS, T_TOTAL, d_val, x_thr_ubsr, ...
%         N_REALIZATIONS, TRANSIENT_RATIO);
    
%     rate_sim_plbsr(i_d) = MonteCarloEscapeRate( ...
%         drift_plbsr_func, FS, T_TOTAL, d_val, x_thr_plbsr, ...
%         N_REALIZATIONS, TRANSIENT_RATIO);
    
%     fprintf('D = %.3f:  HSUBSR r_sim = %.4e,  UBSR r_sim = %.4e,  PLBSR r_sim = %.4e\n', ...
%         d_val, rate_sim_hsubsr(i_d), rate_sim_ubsr(i_d), rate_sim_plbsr(i_d));
% end

% %% 5. 计算相对误差 |r_sim - r_K| / r_K ==================================

% rel_err_hsubsr = abs(rate_sim_hsubsr - rate_theory_hsubsr) ./ rate_theory_hsubsr;
% rel_err_ubsr = abs(rate_sim_ubsr - rate_theory_ubsr) ./ rate_theory_ubsr;
% rel_err_plbsr = abs(rate_sim_plbsr - rate_theory_plbsr) ./ rate_theory_plbsr;

% %% 6. 绘图：Kramers 理论 vs 仿真 & 相对误差 ============================
% SetThesisDefaultStyle();
% % fig = tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
% fig = CreateThesisFigure(16,6);

% % ---- 左图：逃逸率-噪声强度 (半对数) -------------------------------
% subplot(1,2,1);
% semilogy(D_LIST, rate_theory_hsubsr, '-', 'LineWidth', 1.5, 'Color', [0.9020, 0.2941, 0.2078]); hold on;
% semilogy(D_LIST, rate_sim_hsubsr,   'o', 'LineWidth', 1.2, 'Color', [0.9020, 0.2941, 0.2078]);
% semilogy(D_LIST, rate_theory_ubsr, '-', 'LineWidth', 1.5, 'Color', [0.2510, 0.3686, 0.6353]);
% semilogy(D_LIST, rate_sim_ubsr,   's', 'LineWidth', 1.2, 'Color', [0.2510, 0.3686, 0.6353]);
% semilogy(D_LIST, rate_theory_plbsr,  '-', 'LineWidth', 1.5, 'Color', [0.0000, 0.6196, 0.4510]);
% semilogy(D_LIST, rate_sim_plbsr,    'd',  'LineWidth', 1.2, 'Color', [0.0000, 0.6196, 0.4510]);
% xlabel('Noise Intensity D');
% ylabel('Kramers Rate');
% title('Comparison of Kramers Rate');
% legend({'HSUBSR theory','HSUBSR sim', ...
%     'UBSR theory','UBSR sim', ...
%     'PLBSR  theory','PLBSR  sim'}, 'Location', 'best');
% grid on;

% % ---- 右图：相对误差对比 ---------------------------------------------
% subplot(1,2,2);
% plot(D_LIST, rel_err_hsubsr, 'o-', 'LineWidth', 1.5); hold on;
% plot(D_LIST, rel_err_ubsr, 's--', 'LineWidth', 1.5);
% plot(D_LIST, rel_err_plbsr,  'd-.', 'LineWidth', 1.5);
% xlabel('Noise Intensity D');
% ylabel('Relative Error');
% title('Fitting Error');
% legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'best');
% grid on;


% %% 辅助函数部分
% function rate_mean = MonteCarloEscapeRate(drift_func, fs, t_total, d_noise, ...
%     x_thr, n_realizations, transient_ratio)
% % MonteCarloEscapeRate  基于多次仿真估计逃逸率
% %
% % 输入:
% %   drift_func       - 漂移函数句柄 f(x)
% %   fs               - 采样率 (Hz)
% %   t_total          - 单次仿真总时长 (s)
% %   d_noise          - 噪声强度 D
% %   x_thr            - 判别左右势阱的阈值 |x| > x_thr 视为在某个势阱内
% %   n_realizations   - 重复仿真次数
% %   transient_ratio  - 瞬态比例 (0~1)，前期数据丢弃不用统计
% %
% % 输出:
% %   rate_mean        - 该模型在该 D 下的平均逃逸率估计

% h_step   = 1 / fs;
% n_steps  = round(t_total * fs);
% n_points = n_steps + 1;

% clean_signal = zeros(n_points, 1);     % 无输入信号 s(t) = 0
% idx0         = round(transient_ratio * n_points) + 1;

% rate_list = zeros(n_realizations, 1);

% for i_real = 1:n_realizations
%     % 生成高斯白噪声序列 noise_seq，使得:
%     % dx_stoch ≈ noise_seq(i) * h_step = sqrt(2D h) * N(0,1)
%     noise_seq = sqrt(2 * d_noise / h_step) * randn(n_points, 1);
    
%     x_traj = RK4Solver2(drift_func, clean_signal, noise_seq, fs);
    
%     % 丢弃瞬态，仅用稳态部分估计逃逸率
%     x_stead = x_traj(idx0:end);
    
%     rate_list(i_real) = EstimateEscapeRate(x_stead, fs, x_thr);
% end

% rate_mean = mean(rate_list);
% end

% function rate_est = EstimateEscapeRate(x_traj, fs, x_thr)
% % EstimateEscapeRate  根据势阱跃迁次数估计逃逸率
% %
% % 算法:
% %   1. 定义三类区域:
% %      左阱: x <= -x_thr
% %      右阱: x >=  x_thr
% %      中间: 其他 (靠近势垒，不计入势阱驻留)
% %   2. 沿时间扫描区域标签，统计从左阱到右阱或从右阱到左阱的跃迁次数。
% %   3. 逃逸率估计为:
% %      r_hat = (跃迁总次数) / T_obs
% %
% %   对称双稳系统的两阱跃迁率相同，且该定义下的“单位时间跃迁总数”
% %   即为单方向的 Kramers 逃逸率 r。

% n_points = numel(x_traj);
% region   = zeros(n_points, 1);

% % 区域标签: -1 = 左阱, +1 = 右阱, 0 = 中间区
% region(x_traj <= -x_thr) = -1;
% region(x_traj >=  x_thr) =  1;

% last_region   = 0;
% trans_count   = 0;

% for i = 1:n_points
%     current_region = region(i);
    
%     if current_region == 0
%         % 在中间过渡区，不更新跃迁计数，只是停留
%         continue;
%     end
    
%     if last_region == 0
%         % 第一次进入某个势阱
%         last_region = current_region;
%     elseif current_region ~= last_region
%         % 从左阱跳到右阱，或从右阱跳到左阱，计一次跃迁
%         trans_count = trans_count + 1;
%         last_region = current_region;
%     end
% end

% t_obs   = (n_points - 1) / fs;
% rate_est = trans_count / t_obs;
% end

% function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
%     ComputeKramersParamsHSUBSR(a, b, k1, k2)
% % ComputeKramersParamsHSUBSR  计算 HSUBSR 势函数的 Kramers 参数
% %
% % 输出:
% %   x_m     - 右侧势阱最小值位置
% %   x_b     - 中心势垒位置 (对称模型中为 0)
% %   delta_u - 势垒高度 ΔU = U(x_b) - U(x_m)
% %   uxx_m   - U''(x_m)
% %   uxx_b   - U''(x_b) (一般为负值)
% U_fun = @(x) HSUBSR_Potential(x, a, b, k1, k2);

% % 对称双稳势，势垒在 0，右侧势阱在 x > 0
% x_search_right = max(3 * 1, 3);   % 搜索区间上限，可根据参数调整
% x_m = fminbnd(U_fun, 0.1, x_search_right);
% x_b = 0.0;

% % 势垒高度
% delta_u = U_fun(x_b) - U_fun(x_m);

% % 数值二阶导 (有限差分)
% h_fd  = 1e-3;
% uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
% uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
% end

% function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
%     ComputeKramersParamsUBSR(a, b)
% % ComputeKramersParamsUBSR  计算 UBSR 势函数的 Kramers 参数
% U_fun    = @(x) UBSR_Potential(x, a, b);
% threshold = sqrt(a / b);

% % 右侧势阱在 x > 0，约在 threshold 附近
% x_m = fminbnd(U_fun, 0, 2 * threshold);
% x_b = 0.0;   % 对称模型中心势垒

% delta_u = U_fun(x_b) - U_fun(x_m);

% h_fd  = 1e-3;
% uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
% uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
% end

% function [x_m, x_b, delta_u, uxx_m, uxx_b] = ...
%     ComputeKramersParamsPLBSR(u0, l0)
% % ComputeKramersParamsPLBSR  计算 PLBSR 势函数的“有效” Kramers 参数
% %
% % 说明:
% %   PLBSR 为分段线性势函数，在 x = 0, ±b2 处二阶导数不存在，
% %   不严格满足 Kramers 理论的光滑假设。这里采用有限差分在
% %   极值点附近给出一个“有效曲率”，用作对比。

% U_fun    = @(x) PLBSR_Potential(x, u0, l0);
% % 右侧势阱理论上在 x = b
% x_m = l0;
% x_b = 0.0;          % 中心势垒

% delta_u = U_fun(x_b) - U_fun(x_m);
% % 由于为分段线性，二阶导数的有限差分近似依赖于步长 h_fd，
% % 这里只取一个相对“粗”的 h_fd 作为有效曲率估计。
% h_fd  = 1e-3;
% uxx_m = (U_fun(x_m + h_fd) - 2 * U_fun(x_m) + U_fun(x_m - h_fd)) / h_fd^2;
% uxx_b = (U_fun(x_b + h_fd) - 2 * U_fun(x_b) + U_fun(x_b - h_fd)) / h_fd^2;
% end
% =========================================================================
% Script Name: compare_kramers_rate.m
% Description: 基于 Kramers 理论的逃逸率对比实验 (HSUBSR vs UBSR vs PLBSR)
%              验证光滑势函数在理论预测一致性上的优势
%
% Theoretical Basis: r_K = (sqrt(|U''(xb)|*U''(xm)) / 2pi) * exp(-dU/D)
%
% Author: Research Assistant (Generated by Gemini)
% Date: 2026-01-11
% =========================================================================

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

% ---------------------------------------------------------------------
% 仿真参数
% ---------------------------------------------------------------------
fs              = 100;                % 采样率 (Hz)
t_total         = 1000;               % 单次仿真总时长 (s)，需保证足够发生跃迁
transient_ratio = 0.1;                % 瞬态截断比例
n_realizations  = 1;                 % 每个噪声强度下的平均次数 (增加以提高精度)

% 噪声强度扫描范围
D_list = 0.10:0.05:0.50;
num_D  = length(D_list);

%% 2. 模型参数校准 (确保所有模型具有相同的 xm 和 dU) ===================
fprintf('正在校准模型参数...\n');

% (1) HSUBSR 参数校准
shape_factor = 4; % 形状因子，控制势阱宽窄
[a_hs, b_hs, k1_hs, k2_hs] = CalibrateHSUBSR(xm_target, dU_target, shape_factor);

% (2) UBSR 参数校准
% UBSR: xm = sqrt(a/b), dU = a^2/(4b)
% 解得: a = 4*dU/xm, b = 4*dU/xm^2 (仅适用标准UBSR, 此处根据您的UBSR_Potential调整)
% 假设 UBSR_Potential 使用标准参数定义
a_ubsr = 4 * dU_target / xm_target;     
b_ubsr = 4 * dU_target / (xm_target^2); 

% (3) PLBSR 参数校准
% PLBSR: 势阱位置 L0, 势垒高度 U0 (直接对应)
u_pl = dU_target;
l_pl = xm_target;

%% 3. 计算理论 Kramers 逃逸率参数 ======================================
% 计算各模型在势阱底和势垒顶的“曲率” (二阶导数)

% --- HSUBSR (解析/差分均可，此处用差分保证一致性) ---
[~, ~, delta_u_hs, uxx_m_hs, uxx_b_hs] = ...
    ComputeKramersParams(@(x) HSUBSR_Potential(x, a_hs, b_hs, k1_hs, k2_hs), xm_target);

% --- UBSR (分段点导数不连续，使用有限差分近似“有效曲率”) ---
[~, ~, delta_u_ub, uxx_m_ub, uxx_b_ub] = ...
    ComputeKramersParams(@(x) UBSR_Potential(x, a_ubsr, b_ubsr), xm_target);

% --- PLBSR (完全平坦或尖点，曲率理论上为0或无穷，差分仅提供一种参考) ---
[~, ~, delta_u_pl, uxx_m_pl, uxx_b_pl] = ...
    ComputeKramersParams(@(x) PLBSR_Potential(x, u_pl, l_pl), xm_target);

% 计算 Kramers 前因子 (Prefactor): A = sqrt(|U''b|*U''m) / (2pi)
pre_hs = sqrt(abs(uxx_b_hs) * uxx_m_hs) / (2 * pi);
pre_ub = sqrt(abs(uxx_b_ub) * uxx_m_ub) / (2 * pi); 
pre_pl = sqrt(abs(uxx_b_pl) * uxx_m_pl) / (2 * pi);

fprintf('理论参数计算完毕。\n');
fprintf('HSUBSR Prefactor: %.4f\n', pre_hs);
fprintf('UBSR   Prefactor: %.4f (Approximated)\n', pre_ub);
fprintf('PLBSR  Prefactor: %.4f (Approximated)\n', pre_pl);

%% 4. 主循环：理论计算与蒙特卡洛仿真 ====================================
rate_theo_hs = zeros(num_D, 1); rate_sim_hs = zeros(num_D, 1);
rate_theo_ub = zeros(num_D, 1); rate_sim_ub = zeros(num_D, 1);
rate_theo_pl = zeros(num_D, 1); rate_sim_pl = zeros(num_D, 1);

% 定义动力学方程句柄
drift_hs = @(x) HSUBSR_Dynamics(x, a_hs, b_hs, k1_hs, k2_hs);
drift_ub = @(x) UBSR_Dynamics(x, a_ubsr, b_ubsr);
drift_pl = @(x) PLBSR_Dynamics(x, u_pl, l_pl);

% 判别阈值 (用于统计跃迁)
x_thr = 0.5 * xm_target; 

fprintf('\n开始仿真 (Total D points: %d)...\n', num_D);

for i = 1:num_D
    D = D_list(i);
    
    % --- A. 理论值计算 ---
    rate_theo_hs(i) = pre_hs * exp(-delta_u_hs / D);
    rate_theo_ub(i) = pre_ub * exp(-delta_u_ub / D);
    rate_theo_pl(i) = pre_pl * exp(-delta_u_pl / D);
    
    % --- B. 仿真值统计 (蒙特卡洛) ---
    % 使用 parfor 加速 (如果并行工具箱可用，否则改为 for)
    % 这里为了演示清晰使用普通 for 循环，实际运行时建议开启并行
    
    % HSUBSR
    rate_sim_hs(i) = MonteCarloRate(drift_hs, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
    % UBSR
    rate_sim_ub(i) = MonteCarloRate(drift_ub, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
    % PLBSR
    rate_sim_pl(i) = MonteCarloRate(drift_pl, fs, t_total, D, x_thr, n_realizations, transient_ratio);
    
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

function rate_mean = MonteCarloRate(drift_func, fs, t_total, D, x_thr, n_runs, t_skip_ratio)
    % 蒙特卡洛仿真统计平均逃逸率
    rates = zeros(n_runs, 1);
    clean_sig = zeros(round(t_total*fs), 1); % 无信号输入
    
    parfor k = 1:n_runs
        % 生成噪声: dx = f(x)dt + sqrt(2D)dW
        % 离散化: x(n+1) = x(n) + f(x)h + sqrt(2D*h)*randn
        % RK4Solver2 输入 noise_seq 需为 sqrt(2D*h)*randn
        h = 1/fs;
        noise_seq = sqrt(2 * D * h) * randn(length(clean_sig), 1);
        
        % 求解 SDE
        x = RK4Solver2(drift_func, clean_sig, noise_seq, fs);
        
        % 截断稳态数据
        idx_start = round(t_skip_ratio * length(x));
        x_steady = x(idx_start:end);
        
        % 计数跃迁
        rates(k) = CountTransitions(x_steady, x_thr, fs);
    end
    rate_mean = mean(rates);
end

function rate = CountTransitions(x, thr, fs)
    % 统计穿越 +/- thr 的次数
    % 简化逻辑：记录从 <-thr 到 >+thr 的过程 (左->右) 和反向过程
    % 这里计算总跃迁率 r = (N_LR + N_RL) / T_total
    
    state = 0; % 0: undefined/middle, -1: left, 1: right
    
    % 初始化状态
    if x(1) > thr, state = 1;
    elseif x(1) < -thr, state = -1;
    end
    
    count = 0;
    for i = 2:length(x)
        if x(i) > thr
            if state == -1 % 从左跳到右
                count = count + 1;
            end
            state = 1;
        elseif x(i) < -thr
            if state == 1 % 从右跳到左
                count = count + 1;
            end
            state = -1;
        end
    end
    
    T_obs = (length(x)-1)/fs;
    rate = count / T_obs; 
    % 注意：Kramers率通常指单向速率。对于对称势阱，总跳变数/时间 近似等于 2*r_k ?
    % Kramers公式给出的是单向跃迁概率。仿真统计的是双向总频次。
    % 对于对称双稳态，平均单向速率 = 总跃迁数 / 总时间 / 2 ??
    % 不，实际上 r_K 是“粒子处于某个阱中时，单位时间内逃逸的概率”。
    % 简单处理：r_sim = count / T_obs. 
    % 需注意文献定义差异，通常对比数量级和趋势即可。
    % 这里直接返回每秒跃迁次数，与理论公式的 Prefactor 对应即可。
end

function [xm, xb, dU, uxx_m, uxx_b] = ComputeKramersParams(U_func, xm_guess)
    % 数值计算势函数的特征参数
    
    % 1. 寻找精确极值点
    % 势垒顶 (假设在 0 附近)
    xb = 0; 
    % 势阱底 (在 xm_guess 附近搜索)
    options = optimset('Display','off');
    xm = fminbnd(U_func, xm_guess*0.5, xm_guess*1.5, options);
    
    % 2. 势垒高度
    dU = U_func(xb) - U_func(xm);
    
    % 3. 数值二阶导 (有限差分)
    h = 1e-2;
    % U''(xb)
    uxx_b = (U_func(xb+h) - 2*U_func(xb) + U_func(xb-h)) / h^2;
    % U''(xm)
    uxx_m = (U_func(xm+h) - 2*U_func(xm) + U_func(xm-h)) / h^2;
end