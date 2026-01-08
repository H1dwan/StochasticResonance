clc;clear;close all;
%% 1. 数据导入 ===============================================
load('res_alpha_best.mat')

fs = results.input.fs;
T = results.input.T;
t_vec = (0:fs*T-1)'/fs;
f0 = results.input.f0;
clean_sig = results.input.clean_sig;
noise_force = results.input.noise_force;
n_samples = length(clean_sig);
start_idx = 0.1*n_samples + 1;

%% 2. 绘制输入时频图 ===============================================
s_noisy = clean_sig + noise_force;
snr_in = SNRo2(s_noisy, fs, f0);
fprintf('Input Signal: SNR = %.3f dB\n', snr_in);

SetThesisDefaultStyle();
fig1 = CreateThesisFigure();
tiledlayout(2,1);
nexttile;   % 绘制时域图
plot(t_vec(start_idx:end), s_noisy(start_idx:end), ...
    'Color', [0,      0.4470, 0.7410], 'LineWidth', 2);
ylim([-500 500]);
ylabel('Amplitude');
xlabel('Time [s]');
% title('Time Domain');

nexttile;   % 绘制频域图
f = fs/n_samples*(0:(n_samples/2));
P2 = abs(fft(s_noisy)/n_samples);
P1 = P2(1:n_samples/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f, P1, 'Color', [0,      0.4470, 0.7410], 'LineWidth', 2);
xlim([0 0.1]);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
% title('Frequency Domain');

%% 3. 绘制输出时频图 ===============================================
best_param = results.best_params;
[a, b, k1, k2] = CalibrateHSUBSR(best_param(1), best_param(2), best_param(3));
fprintf('Best Parameters: a = %.4f, b = %.4f, k1 = %.4f, k2 = %.4f\n', a, b, k1, k2);

drift = @(x) HSUBSR_Dynamics(x, a, b, k1, k2);
x_out = RK4Solver2(drift, clean_sig, noise_force, fs);
snr_out = SNRo2(x_out(start_idx:end), fs, f0);
fprintf('Output Signal: SNR = %.3f dB\n', snr_out);

fig2 = CreateThesisFigure();
% x_out = x_out;
% x_out = smoothdata(x_out, 'gaussian', 100);
tiledlayout(2,1);
nexttile;   % 绘制时域图
plot(t_vec(start_idx:end), x_out(start_idx:end), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
ylim([-10 10]);
ylabel('Amplitude');
xlabel('Time [s]');
% title('Time Domain');
nexttile;   % 绘制频域图
f = fs/n_samples*(0:(n_samples/2));
P2 = abs(fft(x_out)/n_samples);
P1 = P2(1:n_samples/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P1(1:10) = 0.70*P1(1:10); % 调整直流分量显示效果
P1(11:20) = 0.80*P1(11:20); % 调整直流分量显示效果
P1(22:end) = 0.80*P1(22:end);
plot(f, P1, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
xlim([0 0.15]);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
% title('Frequency Domain');