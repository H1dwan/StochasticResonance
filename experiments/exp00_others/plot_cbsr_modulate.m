%% Demo_TiltedPotential_CBSR.m
% 演示经典双稳系统(CBSR)在周期驱动下的时变倾斜势能展示
% 功能：通过 3D 曲面和 2D 切片展示势能场随时间的动态演化
%
% 参考文献：
% [1] Gammaitoni L, et al. Stochastic Resonance[J]. Reviews of Modern Physics, 1998.

clear; clc; close all;

%% 1. 常量定义
A_PARAM = 1.0;          % 势能参数 a
B_PARAM = 1.0;          % 势能参数 b
SIGNAL_AMP = 0.3;       % 周期信号幅度 A (弱信号)
SIGNAL_FREQ = 0.01;     % 周期信号频率 f (Hz)
X_RANGE = 2.0;          % X轴显示范围 [-X, X]
TIME_PERIODS = 1;       % 展示的周期数
NUM_STEPS = 1000;        % 时间轴采样点数

%% 2. 坐标网格生成
x_vec = linspace(-X_RANGE, X_RANGE, 200);
t_vec = linspace(0, TIME_PERIODS / SIGNAL_FREQ, NUM_STEPS);
[X_mesh, T_mesh] = meshgrid(x_vec, t_vec);

%% 3. 计算等效时变势能 U(x, t)
% 公式：U(x, t) = -0.5*a*x^2 + 0.25*b*x^4 - A*x*cos(2*pi*f*t)
external_signal = SIGNAL_AMP * cos(2 * pi * SIGNAL_FREQ * T_mesh);
potential_v = -0.5 * A_PARAM * X_mesh.^2 + 0.25 * B_PARAM * X_mesh.^4;
U_total = potential_v - X_mesh .* external_signal;

%% 4. 可视化 - 3D 势能曲面演化
% figure('Color', 'w', 'Name', 'CBSR 时变倾斜势能 3D 展示');
% mesh(X_mesh, T_mesh, U_total);
% colormap(jet);
% xlabel('位移 x');
% ylabel('时间 t (s)');
% zlabel('势能 U(x,t)');
% title(['CBSR 时变势能场 (A = ', num2str(SIGNAL_AMP), ')']);
% view(-30, 45);
% grid on;

%% 5. 可视化 - 关键相位下的 2D 切片
% 选取四个关键相位：0, pi/2, pi, 3pi/2
phase_indices = round(linspace(1, NUM_STEPS/TIME_PERIODS, 5));
phase_labels = {'0', '\pi/2', '\pi', '3\pi/2'};
colors = lines(4);



SetThesisDefaultStyle();
color_order = get(0, 'DefaultAxesColorOrder');


for i = 1:4
    idx = phase_indices(i);
    CreateThesisFigure(6, 4.5);
    tiledlayout(1, 1,'Padding','tight','TileSpacing','tight');
    ax = nexttile;
    ax.Position = [0.01 0.01 0.98 0.98];
    plot(x_vec, U_total(idx, :), 'LineWidth', 2);
    xlabel('$x$');
    ylabel('$U(x)$');
    % legend(sprintf('$\phi$ = %s', phase_labels{i}), 'Location', 'northeast');
end

% % 绘制未偏置的原始势能 (A=0)
% v_static = -0.5 * A_PARAM * x_vec.^2 + 0.25 * B_PARAM * x_vec.^4;
% plot(x_vec, v_static, '--k', 'LineWidth', 1.5, 'DisplayName', '原始势能 (A=0)');

% xlabel('位移 x');
% ylabel('势能 U(x)');
% title('周期信号驱动下的势阱倾斜效应');
% legend('Location', 'northeast');
% grid on;
% set(gca, 'FontSize', 12);

%% 6. 动画演示 (可选)
% fprintf('正在生成动态演示...\n');
% figure('Color', 'w', 'Name', '势能场动态演化');
% h_plot = plot(x_vec, U_total(1,:), 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);
% axis([-X_RANGE X_RANGE min(U_total(:))-0.2 max(U_total(:))+0.2]);
% xlabel('位移 x');
% ylabel('等效势能 U(x,t)');
% grid on;

% for k = 1:NUM_STEPS
%     set(h_plot, 'YData', U_total(k, :));
%     title(['时变势能演化 t = ', num2str(t_vec(k), '%.2f'), ' s']);
%     drawnow;
%     pause(0.05);
% end