clc; clear;

load('snr_13_005Hz.mat');
load('mpe_15_005Hz.mat');
load('ami_13_005Hz.mat');

D_list = 0.05:0.01:0.45;

figure;
plot(D_list, snr_mean, 'o-');

figure;
plot(D_list, mpe_mean, 'o-');

figure;
plot(D_list, ami_mean, 'o-');