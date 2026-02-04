% =========================================================================
% Description: 取'res_opt_once_3.mat'中的输入和固定参数
%              取'res_opt_once_2.mat'中的逐层参数
%              取'res_opt_once_best.mat'中的全局参数
%
% Author: LiuShuang
% Created: 2026-02-03
% Last Modified: 2026-02-03
%
% Usage: How to use this script
%
% Input:
%   input_description
% Output:
%   output_description
% =========================================================================

clc; clear; close all;


load('res_opt_once_best.mat');
noise = results.noise;
signal = results.sig;
N = length(noise);
steay_start = round(0.1 * N) + 1;

params_fixed = results.best_params.Fixed;
params_layer = results.best_params.Layer;
params_global = results.best_params.Global;

fs = 5;
f0 = 0.01;

fprintf('输入信号 SNR: %.4f dB\n', SNRo2(noise+signal, fs, f0));
% Plot_Time_Frequency(noise+signal, fs, length(noise), "LineWidth", 1.0);

% csr_fixed = Build_Cascade_SR(3, [params_fixed params_fixed params_fixed]);
% output_fixed = Process_Cascade_SR(noise+signal, csr_fixed, fs);
% fprintf('固定参数级联SR输出信号 SNR: %.4f dB\n', SNRo2(output_fixed(steay_start:end), fs, f0));
% Plot_Time_Frequency(output_fixed, fs, length(output_fixed), "LineWidth", 1.0);

% csr_layer = Build_Cascade_SR(3, params_layer);
% output_layer = Process_Cascade_SR(noise+signal, csr_layer, fs);
% fprintf('分层参数级联SR输出信号 SNR: %.4f dB\n', SNRo2(output_layer(steay_start:end), fs, f0));
% Plot_Time_Frequency(output_layer, fs, length(output_layer), "LineWidth", 1.0);

csr_global = Build_Cascade_SR(3, params_global);
output_global = Process_Cascade_SR(noise+signal, csr_global, fs);
fprintf('全局参数级联SR输出信号 SNR: %.4f dB\n', SNRo2(output_global(steay_start:end), fs, f0));
Plot_Time_Frequency(output_global, fs, length(output_global), "LineWidth", 1.0);



function cascade_system = Build_Cascade_SR(num_layers, params)
% Build_Cascade_SR 构建级联随机共振系统
%
% 输入参数:
%   num_layers      - 级联层数
%   params          - 1x(3*num_layers) 的参数数组
%                     [a1, b1, c1, a2, b2, c2, ..., an, bn, cn]
%                     其中 a, b, c 分别对应各层的三个参数
%
% 输出参数:
%   cascade_system  - 包含所有子系统配置的结构体
% 示例:
%   csr_system = Build_Cascade_SR(3, [a1, b1, c1, a2, b2, c2, a3, b3, c3]);

cascade_system = struct();
cascade_system.num_layers = num_layers;
cascade_system.layers = [];

% 验证参数数组长度
if length(params) ~= 3 * num_layers
    error('参数数组长度错误: 期望 %d 个参数，实际 %d 个', 3*num_layers, length(params));
end

% 为每一层提取参数并构建子系统
for i = 1:num_layers
    layer = struct();
    layer.layer_idx = i;
    
    % 从扁平数组中提取第i层的参数
    idx = (i - 1) * 3;
    [layer.a, layer.b, layer.k1, layer.k2] = CalibrateHSUBSR(...
        params(idx + 1), params(idx + 2), params(idx + 3));
    
    cascade_system.layers = [cascade_system.layers; layer];
end
end

function output = Process_Cascade_SR(input_signal, cascade_system, fs)
% Process_Cascade_SR 通过级联SR系统处理信号

output = input_signal;

for i = 1:cascade_system.num_layers
    layer = cascade_system.layers(i);
    drift_func = @(x) HSUBSR_Dynamics(x, layer.a, layer.b, layer.k1, layer.k2);
    % 调用HSUBSR处理该层
    output = RK4Solver(drift_func, output, 1/fs);
    output = output - mean(output);  % 去直流
end
end