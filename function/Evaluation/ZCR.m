function zcr = ZCR(x, fs)
% 计算原始信号的 ZCR，作为真实的 Kramers 速率统计量。
s = int8(x(:) >= 0);
N = length(x);
num_flips = sum(abs(diff(s)));
duration_sec = N / fs;
zcr = num_flips / duration_sec;
end
