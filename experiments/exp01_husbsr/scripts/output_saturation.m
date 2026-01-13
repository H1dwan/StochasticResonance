% =========================================================================
% Description: 验证 CBSR 在纯周期信号驱动下的输出饱和现象
%
% Author: LiuShuang
% Created: 2026-01-13
% Last Modified: 2026-01-13
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;


%% 1. 系统与仿真参数设置
% -------------------------------------------------------------------------
a = 1;                  % CBSR 参数 a
b = 1;                  % CBSR 参数 b
xm = sqrt(a/b);         % 理论势阱位置 (无信号时)

fs = 10;               % 采样频率 Hz (需满足 fs >> f0)
f0 = 0.01;              % 信号频率 Hz (满足绝热近似条件 f0 << 1)
period = 1/f0;          % 信号周期 s
n_periods = 10;         % 仿真周期数 (保证进入稳态)
T = n_periods * period; % 总时长 s
N = round(fs * T);      % 总点数
t = (0:N-1)'/fs;        % 时间向量

% 定义动力学漂移项句柄
drift = @(x) CBSR_Dynamics(x, a, b);

%% 2. 扫描输入幅值 A，计算输出幅值
% -------------------------------------------------------------------------
A_list = 0.0:0.01:2.0;   % 输入幅值扫描范围 (包含线性区和饱和区)
Out_Amp_list = zeros(size(A_list));

fprintf('开始仿真: 扫描输入幅值 A = %.1f ~ %.1f ...\n', A_list(1), A_list(end));

% 只需要计算一次零噪声 (纯确定性仿真)
noise_seq = zeros(N, 1);

for i = 1:length(A_list)
    Ai = A_list(i);
    
    % 1. 生成输入信号
    s = Ai * sin(2*pi*f0*t);
    
    % 2. 求解系统响应 (RK4)
    x = RK4Solver2(drift, s, noise_seq, fs);
    
    % 3. 截取稳态数据 (取最后 5 个周期)
    steady_points = round(5 * period * fs);
    x_steady = x(end - steady_points + 1 : end);
    
    % 4. 计算输出幅值
    % 方法：计算峰峰值的一半 (Peak-to-Peak / 2)
    % 这比单纯求 max(abs(x)) 更准确，因为信号可能存在直流偏置
    Out_Amp_list(i) = (max(x_steady) - min(x_steady)) / 2;
    
    % 进度条 (可选)
    if mod(i, 5) == 0
        fprintf('  已处理 A = %.2f, 输出幅值 = %.4f\n', Ai, Out_Amp_list(i));
    end
end

fprintf('仿真完成。\n');

%% 3. 结果可视化
% -------------------------------------------------------------------------
SetThesisDefaultStyle();
fig1 = CreateThesisFigure(8, 6, 3);
% 绘制实际响应曲线
plot(A_list, Out_Amp_list, '-x', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
xlabel('$A_{in}$');
ylabel('$A_{out}$');
ylim([0, 2.5]);

% 标记 0.2 和 1.2 两个输入幅值对应的输出响应
plot(0.2, Out_Amp_list(A_list==0.2), 'ko', 'MarkerSize', 8, 'DisplayName', '$A_{in}=0.2$');
plot(1.2, Out_Amp_list(A_list==1.2), 'ko', 'MarkerSize', 8, 'DisplayName', '$A_{in}=1.2$');

fig2 = CreateThesisFigure(4, 3);
s1 = 0.2 * sin(2*pi*f0*t);
x1 = RK4Solver2(drift, s1, noise_seq, fs);
plot(t, x1, 'Color', [0.9020, 0.2941, 0.2078], 'LineWidth', 1.2, 'DisplayName', '$A_{in}=0.2$');
ylim([-2, 2]);
xlabel('$t$');
ylabel('$x(t)$');
legend('Location', 'NorthEast');

fig3 = CreateThesisFigure(4, 3);
s2 = 1.2 * sin(2*pi*f0*t);
x2 = RK4Solver2(drift, s2, noise_seq, fs);
plot(t, x2, 'Color', [0.0000, 0.6196, 0.4510], 'LineWidth', 1.2, 'DisplayName', '$A_{in}=1.2$');
ylim([-2, 2]);
xlabel('$t$');
ylabel('$x(t)$');
legend('Location', 'NorthEast');

%% 4. 辅助分析 (命令行输出)
% -------------------------------------------------------------------------
fprintf('\n=== 结果分析 ===\n');
fprintf('1. 小信号区 (A=%.1f): 输出幅值 = %.4f (近似线性)\n', A_list(1), Out_Amp_list(1));
fprintf('2. 饱和区 (A=%.1f): 输出幅值 = %.4f (增长停滞)\n', A_list(end), Out_Amp_list(end));
fprintf('3. 势阱位置 x_m = %.4f\n', xm);
fprintf('注意观察曲线如何逐渐偏离线性参考线并趋于平缓。\n');