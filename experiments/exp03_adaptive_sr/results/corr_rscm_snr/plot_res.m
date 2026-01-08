clc;clear;close all;

%% 1. 加载数据 ==========================================================
load('res_eta_snr_rscm_corr.mat');
eta_list = results.eta_list;
snr_list = results.snr_list;
J_total_list = results.J_total_list;

%% 2. 计算相关系数矩阵 ==================================================
corr_eta_J      = corr(eta_list, J_total_list);
corr_eta_snr    = corr(eta_list, snr_list);
corr_J_snr      = corr(J_total_list, snr_list);

fprintf('corr(eta, J_total)   = %.3f\n', corr_eta_J);
fprintf('corr(eta, SNR)       = %.3f\n', corr_eta_snr);
fprintf('corr(J_total, SNR)   = %.3f\n', corr_J_snr);

%% 3. 可视化 eta 分布 & eta–J/SNR 散点图 ================================
SetThesisDefaultStyle();
fig = CreateThesisFigure();

subplot(2,2,1);
% histogram(eta_list, 10, 'FaceColor', [0,      0.4470, 0.7410]);
histogram(eta_list, 'BinEdges', 0.2:0.08:1.0, 'FaceColor', [0, 0.4470, 0.7410]);
xticks(0.2:0.2:1.0);
xticklabels({'0.2', '0.4', '0.6', '0.8', '1.0'});
xlabel('$\eta $'); ylabel('counts');
title('$\eta $ distribution');

subplot(2,2,2);
scatter(eta_list, J_total_list, 'filled', 'MarkerFaceColor', [0,      0.4470, 0.7410]);
xlabel('$\eta $'); ylabel('$J_{total}$');
title(sprintf('corr($\\eta $, $J_{total}$) = %.3f', corr_eta_J));

subplot(2,2,3);
scatter(eta_list, snr_list, 'filled', 'MarkerFaceColor', [0,      0.4470, 0.7410]);
xlabel('$\eta $'); ylabel('SNR (dB)');
title(sprintf('corr($\\eta $, SNR) = %.3f', corr_eta_snr));

subplot(2,2,4);
scatter(J_total_list, snr_list, 'filled', 'MarkerFaceColor', [0,      0.4470, 0.7410]);
xlabel('$J_{total}$'); ylabel('SNR (dB)');
title(sprintf('corr($J_{total}$, SNR) = %.3f', corr_J_snr));