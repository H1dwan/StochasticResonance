function [fitness, eta, debug] = RSCMEvaluator(theta_vec, clean_sig, noise_seq, fs, opts)
%RSCMEvaluator 评估 HSUBSR 参数组合的适应度。
%   [fitness, eta, debug] = RSCMEvaluator(theta_vec, clean_sig, noise_seq, fs, opts)
%   在给定信号与噪声输入下，模拟 HSUBSR 动态并提取局部结构指标，
%   通过 AMI 与多尺度排列熵构造综合适应度。
%
% 输入参数
%   theta_vec : [xm, dU, shape_factor]，HSUBSR 结构参数。
%   clean_sig : 原始输入信号列向量。
%   noise_seq : 噪声序列，与 clean_sig 等长。
%   fs        : 采样频率 (Hz)。
%   opts      : 结构体，需包含字段：
%               kappa_eta 共振度 sigmiod 陡峭系数；
%               w_ami, w_mpe 分别为 AMI、MPE 指标的幂次权重。
%
% 输出参数
%   fitness   : 综合适应度（越大越优）。
%   eta       : 基于局部结构的共振度指标，范围 [0,1]。
%   debug     : 调试信息结构体（可选），包含 J_occ/J_ami/J_mpe/eta/scales。

% 解包待优化参数
[a, b, k1, k2] = CalibrateHSUBSR(theta_vec(1), theta_vec(2), theta_vec(3));
% a = theta_vec(1);
% b = theta_vec(2);
% k1 = theta_vec(3);
% k2 = theta_vec(4);

% 1. 仅截取稳态段用于指标计算，剔除前 10% 过渡期
n_samples    = length(clean_sig);
steady_start = round(0.1 * n_samples);

% 2. 构建漂移函数并数值积分得到输出轨迹
drift_func = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x_out = RK4Solver(drift_func, clean_sig + noise_seq, 1/fs);
x_steady = x_out(steady_start+1:end);

% 3. 局部结构特征：阱占据度 + AMI 归一化 + 共振度 eta
[J_occ, J_ami, eta] = LocalStructureFeatures(x_steady, fs, opts.kappa_eta);

% 4. 多尺度排列熵：取归一化序列标准差的反向指标
[scales, ~] = SelectMpeScales(x_steady, fs, 4);
[~, mpnorm, ~, ~] = MultiScalePermEn(x_steady, scales);
% mwpe_index = std(mpnorm(~isnan(mpnorm)));
mwpe_index = mean(mpnorm);
if isnan(mwpe_index)
    J_mpe = 0;
else
    J_mpe = max(0, min(1, 1 - mwpe_index));
end

% 5. 适应度：AMI^w_ami * MPE^w_mpe，均为 [0,1] 内的非负指标
% fitness = J_ami.^opts.w_ami * J_mpe.^opts.w_mpe;
% fitness = (1-J_ami).^opts.w_ami * mwpe_index.^opts.w_mpe;
fitness = -J_ami * J_mpe;  % 取负号以实现最大化

% 调试信息（可选输出）
if nargout >= 3
    debug.J_occ = J_occ;
    debug.J_ami = J_ami;
    debug.J_mpe = J_mpe;
    debug.eta   = eta;
    debug.scales = scales;
end

end

%% 内部辅助函数
function [J_occ, J_ami, eta] = LocalStructureFeatures(x, fs, kappa_eta)
%LocalStructureFeatures 计算局部结构指标：阱占据度、AMI、共振度。
%   输入 x 为单通道序列，fs 为采样率，kappa_eta 为 sigmiod 陡峭系数。

x = x(:);
n_len = length(x);

% 1) 滞回二值化：根据均值±0.3σ 形成两态序列 s(n)
mu_x    = mean(x);
sigma_x = std(x);

th_high = mu_x + 0.3 * sigma_x;
th_low  = mu_x - 0.3 * sigma_x;

s_state = zeros(n_len, 1);
state   = 1;
if x(1) < mu_x
    state = 0;
end

for idx = 1:n_len
    if x(idx) > th_high
        state = 1;
    elseif x(idx) < th_low
        state = 0;
    end
    s_state(idx) = state;
end

% 2) 阱占据平衡度 J_occ：p1=0.5 时取 1，失衡时下降至 0
p1    = mean(s_state);                 % 右阱占据概率
J_occ = 1 - abs(p1 - 0.5) / 0.5;
J_occ = max(0, min(1, J_occ));

% 3) AMI 指标（调用 AMI2），并做归一化 ami/(ami+c0)
[J_ami, ~] = AMI2(x, fs);
c0         = 0.1;                      % 防止小值放大，控制饱和速率
J_ami_norm = J_ami / (J_ami + c0);
J_ami_norm = max(0, min(1, J_ami_norm));

% 4) 共振度 eta：将 J_occ 与 J_ami_norm 线性融合后经 sigmiod 映射
beta_occ   = 1.0;
beta_ami   = 1.0;

R  = beta_occ * J_occ + beta_ami * J_ami_norm;
R0 = 1.0;                             % sigmiod 中心，可根据实验调整

eta = 1 ./ (1 + exp(-kappa_eta * (R - R0)));
eta = max(0, min(1, eta));

end