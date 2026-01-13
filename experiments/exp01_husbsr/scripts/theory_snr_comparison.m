% ===============================================================================
% Description: 对比不同势垒高度、不同势阱位置下 HSUBSR 和 CBSR 的理论信噪比 (SNR) 曲线
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
% ================================================================================

clc; clear; close all;
SetThesisDefaultStyle();

%% 1. 全局参数设置 (Global Parameters)
% 信号参数
A = 0.1;           % 信号幅值 (弱信号)
f = 0.01;          % 信号频率 (Hz, 绝热条件 f << r_k)

% 噪声扫描范围
D_list = 0.0:0.01:1.5;

% 绘图颜色
colors = [ ...
    0,      0.4470, 0.7410; % 深蓝
    0.9020, 0.2941, 0.2078; % 朱红
    0.0000, 0.6196, 0.4510; % 墨绿
    ];

%% 2. 实验一：固定 xm，改变 dU 比较 HSUBSR 和 CBSR 的 SNR 曲线
CreateThesisFigure();
layout = tiledlayout(1,1);       %分区作图
layout.Padding = 'compact';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

dU_list = [0.25, 0.50, 0.75];       % 不同势垒高度
xm_fixed = 1.0;                  % 固定势阱位置
shape = 50;

legend_str = cell(1, 2 * length(dU_list));

for i = 1:length(dU_list)
    dU = dU_list(i);
    
    % --- Model A: HSUBSR 理论计算 ---
    % 1. 参数校准
    [a_h, b_h, k1, k2] = CalibrateHSUBSR(xm_fixed, dU, shape);
    
    % 2. 计算特征频率
    wb_h = sqrt(a_h*k1 - b_h*k2); % 势垒顶曲率频率
    term_w0 = b_h*k2*(sech(k2*xm_fixed))^2 - a_h*k1*(sech(k1*xm_fixed))^2;
    w0_h = sqrt(abs(term_w0));    % 势阱底曲率频率
    
    % 3. 计算 SNR 曲线
    snr_hsubsr = calculate_snr(D_list, w0_h, wb_h, dU, xm_fixed, A, f);
    
    % --- Model B: CBSR 理论计算 ---
    % 1. 参数校准 (a, b)
    % xm = sqrt(a/b), dU = a^2/(4b) => a = 4*dU/xm^2
    a_c = 4 * dU / (xm_fixed^2);
    % b_c = 4 * dU / (xm_fixed^4); (计算r_k不需要b，只需要a)
    
    % 2. 计算特征频率
    % U''(0) = -a => wb = sqrt(a)
    % U''(xm) = 2a => w0 = sqrt(2a)
    wb_c = sqrt(a_c);
    w0_c = sqrt(2 * a_c);
    
    % 3. 计算 SNR 曲线
    snr_cbsr = calculate_snr(D_list, w0_c, wb_c, dU, xm_fixed, A, f);
    
    color = colors(i,:);
    % --- 绘图 ---
    plot(D_list, snr_hsubsr, '-' , 'LineWidth', 2, 'Color', color); hold on;
    plot(D_list, snr_cbsr, '--', 'LineWidth', 2, 'Color', color);
    
    legend_str{2*i-1} = sprintf('HSUBSR ($\\Delta U=%.2f$)', dU);
    legend_str{2*i} = sprintf('CBSR ($\\Delta U=%.2f$)', dU);
end

% 图形美化
xlabel('$D$');
ylabel('$\mathrm{SNR}$');
legend(legend_str, 'Location', 'Northeast');

%% 3. 实验二：固定 dU，改变 xm 比较 HSUBSR 和 CBSR 的 SNR 曲线
CreateThesisFigure();
layout = tiledlayout(1,1);       %分区作图
layout.Padding = 'compact';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

xm_list = [1.00, 2.00, 3.00];       % 不同势阱位置
dU_fixed = 0.25;                    % 固定势垒高度
shape = 50;

legend_str = cell(1, 2 * length(xm_list));
for i = 1:length(xm_list)
    xm = xm_list(i);
    
    % --- Model A: HSUBSR 理论计算 ---
    % 1. 参数校准
    [a_h, b_h, k1, k2] = CalibrateHSUBSR(xm, dU_fixed, shape);
    
    % 2. 计算特征频率
    wb_h = sqrt(a_h*k1 - b_h*k2); % 势垒顶曲率频率
    term_w0 = b_h*k2*(sech(k2*xm))^2 - a_h*k1*(sech(k1*xm))^2;
    w0_h = sqrt(abs(term_w0));    % 势阱底曲率频率
    
    % 3. 计算 SNR 曲线
    snr_hsubsr = calculate_snr(D_list, w0_h, wb_h, dU_fixed, xm, A, f);
    
    % --- Model B: CBSR 理论计算 ---
    % 1. 参数校准 (a, b)
    % xm = sqrt(a/b), dU = a^2/(4b) => a = 4*dU/xm^2
    a_c = 4 * dU_fixed / (xm^2);
    % b_c = 4 * dU / (xm_fixed^4); (计算r_k不需要b，只需要a)
    
    % 2. 计算特征频率
    % U''(0) = -a => wb = sqrt(a)
    % U''(xm) = 2a => w0 = sqrt(2a)
    wb_c = sqrt(a_c);
    w0_c = sqrt(2 * a_c);
    
    % 3. 计算 SNR 曲线
    snr_cbsr = calculate_snr(D_list, w0_c, wb_c, dU_fixed, xm, A, f);
    
    color = colors(i,:);
    % --- 绘图 ---
    plot(D_list, snr_hsubsr, '-' , 'LineWidth', 2, 'Color', color); hold on;
    plot(D_list, snr_cbsr, '--', 'LineWidth', 2, 'Color', color);
    
    legend_str{2*i-1} = sprintf('HSUBSR ($x_m=%.2f$)', xm);
    legend_str{2*i} = sprintf('CBSR ($x_m=%.2f$)', xm);
end

% 图形美化
xlabel('$D$');
ylabel('$\mathrm{SNR}$');
legend(legend_str, 'Location', 'Northeast');

%% 4. 辅助函数 (Local Functions)
% -----------------------------------------------------------

function snr = calculate_snr(D_vec, w0, wb, dU, xm, A, ~)
% 通用 SNR 计算函数 (Based on McNamara & Wiesenfeld / Gammaitoni)
snr = zeros(size(D_vec));
for k = 1:length(D_vec)
    D = D_vec(k);
    % Kramers Rate
    r_k = (w0 * wb / (2*pi)) * exp(-dU / D);
    
    % Adiabatic SNR part
    snr_adiabatic = (xm^2 * A^2 * r_k) / (4 * D^2);
    
    % Frequency Correction (Lorentzian Roll-off)
    % correction = 1 / (1 + (pi * f / r_k)^2);
    
    correction = 1;
    
    snr(k) = snr_adiabatic * correction;
end
end