clc; clear; close all;

%% 1. 加载数据 =========================================================================
load('res_white_noise_comp.mat'); 

fs = results.input.fs;
f0 = results.input.f0;
x_ubsr = results.outputs.UBSR;
x_plbsr = results.outputs.PLBSR;
x_hs = results.outputs.HSUBSR;

%% 2. 打印输出 =========================================================================
fprintf('White Noise Performance Results:\n');
fprintf('UBSR: %.4f\n', SNRo2(x_ubsr, fs, f0));
fprintf('PLBSR: %.4f\n', SNRo2(x_plbsr, fs, f0));
fprintf('HSUBSR: %.4f\n', SNRo2(x_hs, fs, f0));