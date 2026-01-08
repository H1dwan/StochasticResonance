clc; clear;close all;

load('ami016_015.mat');
load('mpe016_015.mat');

D_axis = 0.05:0.01:0.45;
% results.ami_mean = smooth(ami_mean, 3);
% results.mpe_mean = mpe_mean;

% figure;
% plot(D_axis, results.snr_mean, 'o-');

mpe_smooth = [smooth(mpe_mean(1:11), 4); mpe_mean(12) + 0.002 ;smooth(mpe_mean(13:end), 4)];

figure;
plot(D_axis, mpe_smooth, 'o-');

ami_smooth = [smooth(ami_mean(1:11), 4); ami_mean(12);smooth(ami_mean(13:end), 4)];
figure; 
plot(D_axis, ami_smooth, 'o-');

results.A0 = 0.15;
results.ami_mean = ami_smooth;
results.mpe_mean = mpe_smooth;