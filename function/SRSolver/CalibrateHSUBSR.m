function [a, b, k1, k2] = CalibrateHSUBSR(xm, dU, shape)
k2 = 1.0 / xm;          % 尾部衰减尺度
k1 = shape * k2;        % 核心锐度尺度 (越大越接近 UBSR)

% 求解强度参数 a, b 以匹配 dU 和 xm
% 方程组:
% 1. f(xm) = 0  -> a*tanh(k1*xm) = b*tanh(k2*xm) -> a = b * R
% 2. U(0)-U(xm) = dU

R = tanh(k2*xm) / tanh(k1*xm);

% dU = - [ (b/k2)*ln(cosh(k2*xm)) - (a/k1)*ln(cosh(k1*xm)) ]
%    = b * [ (R/k1)*ln(cosh(k1*xm)) - (1/k2)*ln(cosh(k2*xm)) ]

term = (R/k1)*log(cosh(k1*xm)) - (1/k2)*log(cosh(k2*xm));

if term <= 0
    error('Shape factor too small, cannot form bistable potential');
end

b = dU / term;
a = b * R;
end