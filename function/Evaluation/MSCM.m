function [J_mscm, J_mwpe, J_ami] = MSCM(x, fs)
% MSCM 计算改进统计复杂度 (MSCM)、多尺度加权置换熵 (MWPE) 和自互信息 (AMI)
out = AdaptiveScalesACF(x, fs);
[~, mpnorm, ~, ~] = MultiScalePermEn(x, out.S);
J_mwpe = mean(mpnorm(~isnan(mpnorm)));

[J_ami, ~] = AMI2(x, fs);

J_mscm = J_ami * (1 - J_mwpe);
end