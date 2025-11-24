function x_out = RK4Solver2(drift_func, clean_signal, noise_seq, fs, varargin)
% RK4Solver2 随机共振分步求解器 (Split-Step Method)
%
% [方法论]
% 该求解器采用算子分裂法 (Operator Splitting) 处理随机微分方程 (SDE):
% 1. 确定性部分 (漂移项 + 周期信号): 使用高精度 RK4 积分
% 2. 随机性部分 (噪声项): 使用 Euler-Maruyama 积分 (直接叠加)
% 公式: x_{n+1} = x_n + RK4_Step(Drift + Signal) + Noise_n * h
%
% [优势]
% 相比直接将噪声带入RK4，该方法在数学上更严谨，且保留了RK4处理非线性势阱的稳定性。
%
% [输入参数]
%   drift_func   - 函数句柄，系统固有漂移 f(x) = -U'(x) (不含信号和噪声!)
%   clean_signal - 纯输入信号向量 (Deterministic Signal)
%   noise_seq    - 噪声序列向量 (Stochastic Force)，需与信号等长
%                  注意：这里假设输入的是噪声力 xi(t)，积分后为 xi(t)*h
%   fs           - 采样频率 (Hz)
%   varargin     - 可选: 'InitialCondition' (默认0)
%
% [输出]
%   x_out        - 系统输出轨迹

% --- 1. 参数解析 ---
p = inputParser;
addOptional(p, 'InitialCondition', 0);
parse(p, varargin{:});
x0 = p.Results.InitialCondition;

h = 1/fs; % 时间步长
N = length(clean_signal);

% 确保信号和噪声长度一致
if length(noise_seq) ~= N
    error('SR:DimError', '输入信号和噪声序列长度必须一致');
end

% --- 2. 内存预分配 ---
x_out = zeros(N, 1);
x_out(1) = x0;

% 预计算常数，提升循环效率
h_half = 0.5 * h;
h_sixth = h / 6;

% --- 3. 核心求解循环 (Split-Step) ---
current_x = x0;

for i = 1 : N-1
    % === 步骤 A: 准备数据 ===
    % 取当前时刻的确定性信号值
    % (SR中信号变化通常慢于采样，取 s(i) 或 s(i+0.5) 差异不大，这里取左端点)
    s_val = clean_signal(i);
    s_val_mid = 0.5 * (clean_signal(i) + clean_signal(i+1)); % 中点信号估计(可选，提高精度)
    s_val_next = clean_signal(i+1);
    
    % === 步骤 B: 确定性部分的 RK4 演化 (Deterministic Evolution) ===
    % 计算仅包含 Potential + Signal 的增量
    
    % k1: 基于当前点
    k1 = drift_func(current_x) + s_val;
    
    % k2: 基于中点 (预测)
    x_k2 = current_x + k1 * h_half;
    k2 = drift_func(x_k2) + s_val_mid;
    
    % k3: 基于中点 (修正)
    x_k3 = current_x + k2 * h_half;
    k3 = drift_func(x_k3) + s_val_mid;
    
    % k4: 基于终点
    x_k4 = current_x + k3 * h;
    k4 = drift_func(x_k4) + s_val_next;
    
    % 计算确定性增量 dx_det
    dx_det = (k1 + 2*k2 + 2*k3 + k4) * h_sixth;
    
    % === 步骤 C: 随机性部分的 Euler 叠加 (Stochastic Addition) ===
    % 噪声项直接叠加：积分 int(xi(t) dt) approx xi[i] * h
    % 注意：你的 noise_seq 应该是已经生成好的随机力序列
    dx_stoch = noise_seq(i) * h;
    
    % === 步骤 D: 更新状态 ===
    current_x = current_x + dx_det + dx_stoch;
    x_out(i+1) = current_x;
end
end