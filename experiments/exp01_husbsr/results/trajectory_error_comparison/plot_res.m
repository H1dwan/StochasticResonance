clc; clear; close all;
%% 1. 导入数据 =====================================================
load('rms_err_hsubsr.mat'); % 加载 HSUBSR 误差数据
load('rms_err_ubsr.mat');   % 加载 UBSR 误差
load('rms_err_plbsr.mat');  % 加载 PLBSR 误差
fs_list      = [10, 20, 25, 40, 50, 100, 125, 200]; 

%% 2. 绘图 =========================================================
SetThesisDefaultStyle();
CreateThesisFigure();
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
nexttile;
plot(fs_list, rms_err_hsubsr, 'o--', 'LineWidth', 2); hold on;
plot(fs_list, rms_err_ubsr, 's-.', 'LineWidth', 2);
plot(fs_list, rms_err_plbsr,  'd:', 'LineWidth', 2);
xlabel('$f_s$', 'Interpreter', 'latex');
ylabel('RMS Error');
legend({'HSUBSR','UBSR','PLBSR'}, 'Location', 'northeast');
