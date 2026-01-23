function Plot_Time_Frequency(x, fs, N, options)
% 绘制实信号的单边幅频图
%   该函数用于同时绘制信号的时域波形和频域幅频特性曲线
%
% 输入参数:
%   x - 输入实信号，一维数组
%   fs - 采样频率，单位Hz
%   N - 信号长度，即信号点数
%   options - 绘图选项结构体
%       .LineWidth - 线条宽度，默认为1

% 参数解析
arguments
    x   double
    fs  double
    N   double
    options.LineWidth (1,1) {mustBeNumeric} = 1
end

% 时间轴和频率轴计算
t = (0:N-1) / fs;
f = fs/N*(0:(N/2)); % 频率范围（包含奈奎斯特频率）

% FFT计算和幅值归一化
P2 = abs(fft(x)/N); % 计算fft并归一化纵轴（双边）
P1 = P2(1:N/2+1);   % 截取单边谱
P1(2:end-1) = 2*P1(2:end-1);    % 能量还原（除0Hz和奈奎斯特频率）

% 绘制时域和频域图形
SetThesisDefaultStyle();
CreateThesisFigure(8, 3);
layout = tiledlayout(1,2);   %分区作图
layout.Padding = 'tight';      % 紧凑内边距
layout.TileSpacing = 'tight';    % 紧密的图块间距

nexttile
plot(t, x, 'LineWidth', options.LineWidth);
% ylim([-20 20])
xlabel('Time[s]')
ylabel('Amplitude')
% title('Time Domain')
% set(gca,'FontSize',14,'FontName','Times New Roman');
% xlim([0 N])
xticks(0:200:1000);

nexttile;
plot(f, P1, 'LineWidth', options.LineWidth);
xlim([0 0.3])
xlabel('Frequency[Hz]')
ylabel('Amplitude')
% title('Frequency Domain')
% set(gca,'FontSize',14,'FontName','Times New Roman');

end