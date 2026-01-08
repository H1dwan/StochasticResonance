% =========================================================================
% Description: 验证 HSUBSR 势模型在α稳定噪声下的表现，传统的 CBSR 会因为α噪声的
%              脉冲特性，导致数值越界，而 HSUBSR 凭借线性势阱壁则能较好地抑制
%
% Author: LiuShuang
% Created: 2026-01-07
% Last Modified: 2026-01-07
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;

%% 1. 公共仿真参数 =========================================================
fs         = 5;               % 采样频率 Hz
T          = 2000;            % 信号时长 s
n_samples  = fs * T;
t_vec      = (0:n_samples-1)' / fs;
h  = 1/fs;                    % 步长

% 弱信号参数
A0         = 0.35;             % 弱信号幅值
f0         = 0.01;             % 真实信号频率（仅用于验证）
clean_sig  = A0 * sin(2*pi*f0*t_vec);

% Alpha 稳定噪声参数
alpha_noise = 1.2;      % 特征指数 (1.0 表示显著脉冲)
beta_noise  = 0;        % 对称分布
sigma_noise = 2;        % 噪声强度
mu_noise    = 0;        % 零均值

% 生成 Alpha 稳定噪声
sigma_discrete = sigma_noise * h^(1/alpha_noise);
noise_force = GenerateAlphaStableNoise(alpha_noise, beta_noise, sigma_discrete, mu_noise, n_samples, 1) * 1;

% 混合信号
s_noisy = clean_sig + noise_force; % 这里的 s_noisy 仅用于观察，实际求解器输入分开传

% 绘制时频图
Plot_Time_Frequency(s_noisy, fs, n_samples);
xlim([0 0.2]);

% 计算输入 SNR (参考)
snr_in = SNRo2(s_noisy, fs, f0);
fprintf('输入信号 SNR = %.2f dB \n', snr_in);


%% 2. 优化参数设置 (PSO) =====================================================
fprintf('2. 配置 PSO 优化算法...\n');

% HSUBSR 势结构参数搜索范围，使用 "物理映射搜索策略" (xm, dU, shape_factor)
% input_rms = std(clean_sig+noise_force);
% theta_min = [0.5, 0.5, 1];  % [xm, dU, shape_factor]
% theta_max = [3, 5, 50];

% UBSR 势参数搜索范围
theta_min = [0.1, 0.1];  % [a, b]
theta_max = [10, 10];

% PSO 参数
SearchAgents_no = 20;   % 种群规模
Max_iter = 30;          % 最大迭代次数 (演示用30，实际建议50+)

evaluator = @(x) SNREvaluator2(x, clean_sig, noise_force, fs, f0);
% evaluator = @(x) SNREvaluator(x, clean_sig, noise_force, fs, f0);

%% 3. 执行优化 ===============================================
fprintf('3. 开始 PSO 寻优 ...\n');
tic;
[best_score, best_pos, convergence_curve] = PSO(SearchAgents_no, Max_iter, theta_min, theta_max, ...
    length(theta_min), evaluator);
elapsed_time = toc;

best_snr = -best_score; % 还原为正的 SNR
fprintf('   优化完成！耗时: %.2f 秒\n', elapsed_time);
fprintf('   最优 SNR: %.4f dB\n', best_snr);
fprintf('   最优参数: \n');
fprintf('     param1 = %.4f\n', best_pos(1));
fprintf('     param2 = %.4f\n', best_pos(2));
% fprintf('     shape_factor         = %.4f\n', best_pos(3));

%% 4. 结果验证与绘图 =========================================
fprintf('4. 绘制结果...\n');

% 使用最优参数重新运行一次仿真
opt_a = best_pos(1);
opt_b = best_pos(2);
% opt_xm = best_pos(1);
% opt_dU = best_pos(2);
% opt_shape_factor = best_pos(3);

% [a, b, k1, k2] = CalibrateHSUBSR(opt_xm, opt_dU, opt_shape_factor);
% drift_func = @(x) HESUBSR_Dynamics(x, a, b, k1, k2);
drift_func = @(x) UBSR_Dynamics(x, opt_a, opt_b);
x_opt = RK4Solver2(drift_func, clean_sig, noise_force, fs);

% --- 绘图 ---
figure('Position', [100, 100, 1200, 800], 'Color', 'white');

% 1. 收敛曲线
plot(-convergence_curve, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('Iteration'); ylabel('Max SNR (dB)');
title('PSO Convergence Curve');
grid on;

% 2. 时频域波形对比
% Plot_Time_Frequency(x_opt, fs, length(x_opt));
% xlim([0 0.05]);
% figure;
Plot_Time_Frequency(x_opt, fs, length(x_opt));
xlim([0 0.1]);

%% 5. 保存结果 ==========================================================
results.input.A0 = A0;
results.input.f0 = f0;
results.input.fs = fs;
results.input.T = T;
results.input.clean_sig = clean_sig;
results.input.noise_force = noise_force;
results.curve.pso = convergence_curve;
results.best_params = best_pos;

%% 辅助函数
function fitness = SNREvaluator2(theta, clean_sig, noise_seq, fs, f0)
a = theta(1);
b = theta(2);
drift_func = @(x) UBSR_Dynamics(x, a, b);
x_out = RK4Solver2(drift_func, clean_sig , noise_seq, fs);
snr_val = SNRo2(x_out, fs, f0);
fitness = -snr_val;
end

function x = GenerateAlphaStableNoise(alpha, beta, sigma, mu, m, n)
% GenerateAlphaStableNoise 生成 Alpha 稳定分布随机数
% 基于 Chambers-Mallows-Stuck (CMS) 算法
%
% 用法:
%   x = GenerateAlphaStableNoise(alpha, beta, sigma, mu, m, n)
%
% 输入参数:
%   alpha - 特征指数 (0 < alpha <= 2)
%           描述分布的拖尾厚度（脉冲性）。
%           alpha=2 为高斯分布，alpha<2 表现出脉冲特性。
%   beta  - 偏度参数 (-1 <= beta <= 1)
%           beta=0 为对称分布 (S alpha S)，通常用于随机共振研究。
%   sigma - 尺度参数 (sigma > 0)
%           控制噪声的强度/离散程度。
%   mu    - 位置参数 (实数)
%           分布的中心位置（如中位数），通常设为 0。
%   m, n  - 输出矩阵的维度 (行数 m, 列数 n)
%           如果省略，默认生成一个标量。
%
% 输出参数:
%   x     - 大小为 m x n 的随机数矩阵
%
% 参考文献:
%   [1] Chambers, J. M., Mallows, C. L., & Stuck, B. W. (1976).
%       A method for simulating stable random variables.
%       Journal of the American Statistical Association.
%   [2] Weron, R. (1996). On the stable motion of stock prices.

% --- 1. 参数校验与默认值处理 ---
if nargin < 6, n = 1; end
if nargin < 5, m = 1; end
if nargin < 4, mu = 0; end
if nargin < 3, sigma = 1; end

% 检查 alpha 范围
if alpha <= 0 || alpha > 2
    error('特征指数 alpha 必须在 (0, 2] 范围内');
end
% 检查 beta 范围
if abs(beta) > 1
    error('偏度参数 beta 必须在 [-1, 1] 范围内');
end
% 检查 sigma 范围
if sigma <= 0
    error('尺度参数 sigma 必须大于 0');
end

% --- 2. 生成辅助随机变量 ---
% V: 均匀分布于 (-pi/2, pi/2)
V = (rand(m, n) - 0.5) * pi;
% W: 指数分布，均值为 1
W = exprnd(1, m, n);

% --- 3. 核心算法 (CMS Method) ---

% == 情况 A: alpha = 2 (高斯分布) ==
if alpha == 2
    % 直接使用 Box-Muller 变换或 MATLAB 内置 randn
    % S_2(sigma, 0, mu) 等价于 N(mu, 2*sigma^2)
    % 注意标准定义下: sigma_gauss = sigma * sqrt(2)
    x = sigma * sqrt(2) * randn(m, n) + mu;
    return;
end

% == 情况 B: alpha = 1 (柯西分布类型的特殊情况) ==
if alpha == 1
    term1 = (pi/2 + beta * V) .* tan(V);
    
    % 避免 log(0) 或除以 0 的数值保护
    arg_log = (pi/2 * W .* cos(V)) ./ (pi/2 + beta * V);
    arg_log(arg_log <= 0) = eps; % 防止数值错误
    
    term2 = beta * log(arg_log);
    
    % 标准化随机变量 X ~ S_1(1, beta, 0)
    X = 2/pi * (term1 - term2);
    
    % 尺度变换: Y = sigma * X + 2/pi * beta * sigma * log(sigma) + mu
    % 注意：alpha=1 时的移位项比较特殊
    x = sigma * X + (2/pi * beta * sigma * log(sigma)) + mu;
    return;
end

% == 情况 C: alpha != 1 且 alpha != 2 (一般情况) ==
% 计算常量 B_alpha_beta
B = atan(beta * tan(pi * alpha / 2)) / alpha;

% 计算常量 S_alpha_beta
S = (1 + (beta * tan(pi * alpha / 2))^2)^(1 / (2 * alpha));

% 分步计算公式
num = sin(alpha * (V + B));               % 分子
den = cos(V).^(1 / alpha);                % 分母部分 1

% 括号内的项
term_inner = cos(V - alpha * (V + B)) ./ W;
term_inner(term_inner <= 0) = eps;        % 数值保护

term2 = term_inner.^((1 - alpha) / alpha); % 第二部分

% 标准化随机变量 Z ~ S_alpha(1, beta, 0)
Z = S * (num ./ den) .* term2;

% --- 4. 尺度与位置变换 ---
% x ~ S_alpha(sigma, beta, mu)
x = sigma * Z + mu;

end