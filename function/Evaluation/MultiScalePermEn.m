function [mpe, mpnorm, mcpe, details] = MultiScalePermEn(sig, scales, varargin)
% MultiScalePermEn  多尺度排列熵 (Multiscale Permutation Entropy)
%
% 用法:
%   [mpe, mpnorm, mcpe] = MultiScalePermEn(sig, scales, 'Name', Value, ...)
%
% 必选参数:
%   sig     - 输入信号 (列/行向量)
%   scales  - 尺度向量，例如 1:20 或 [1 2 4 8]
%
% 可选参数 (传递给 PermEn):
%   'm'        - 嵌入维度，默认 3
%   'tau'      - 时间延迟，默认 1
%   'Logx'     - 对数底，0=e，默认 2
%   'Norm'     - 归一化开关，默认 false
%   'Typex'    - 变体类型 (见 PermEn)，默认 'none'
%   'tpx'      - 变体参数 (可留空)
%
% 输出:
%   mpe     - 每个尺度的排列熵 (取 PermEn 返回的最高维度项)
%   mpnorm  - 归一化排列熵
%   mcpe    - 条件排列熵 (对应最高维度)
%   details - 结构体，包含每个尺度完整的 Perm/Pnorm/cPE 向量
%
% 依赖: 需要 EntropyHub 的 PermEn 函数已在路径中。

if nargin < 2 || isempty(scales)
    scales = 1:10;
end

p = inputParser;
p.addParameter('m', 3, @(x) isnumeric(x) && isscalar(x) && x > 1);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('Logx', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Norm', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Typex', 'weighted', @(x) ischar(x) || isstring(x));
p.addParameter('tpx', [], @(x) isscalar(x) && (x > 0));
p.parse(varargin{:});
args = p.Results;

sig = sig(:); % 转为列向量

% 确保尺度为正整数向量，保持用户提供的顺序
if ~isvector(scales) || any(scales <= 0) || any(scales ~= floor(scales))
    error('scales 必须为正整数向量，例如 1:20 或 [1 2 4]');
end
scales = scales(:); % 列向量以便索引
S = numel(scales);

mpe    = nan(S, 1);
mpnorm = nan(S, 1);
mcpe   = nan(S, 1);

% 记录每个尺度完整结果以便调试/分析
perm_all  = cell(S, 1);
pnorm_all = cell(S, 1);
cpe_all   = cell(S, 1);

for idx = 1:S
    s = scales(idx);
    cg = coarse_grain(sig, s);
    if numel(cg) < (args.m - 1) * args.tau + 2
        % 数据不足以构造嵌入，保持 NaN
        fprintf('尺度 %d 下数据点不足，跳过计算。\n', s);
        continue;
    end
    
    % 调用 PermEn。若 tpx 为空则不传递该参数，避免覆盖默认。
    if isempty(args.tpx)
        [perm, pnorm, cpe] = PermEn(cg, 'm', args.m, 'tau', args.tau, ...
            'Logx', args.Logx, 'Norm', args.Norm, 'Typex', args.Typex);
    else
        [perm, pnorm, cpe] = PermEn(cg, 'm', args.m, 'tau', args.tau, ...
            'Logx', args.Logx, 'Norm', args.Norm, 'Typex', args.Typex, 'tpx', args.tpx);
    end
    
    % 取最高嵌入维度的值作为该尺度的指标
    mpe(idx)    = perm(end);
    mpnorm(idx) = pnorm(end);
    if ~isempty(cpe)
        mcpe(idx) = cpe(end);
    end
    
    perm_all{idx}  = perm;
    pnorm_all{idx} = pnorm;
    cpe_all{idx}   = cpe;
end

details = struct('scales', scales, 'perm', {perm_all}, 'pnorm', {pnorm_all}, 'cpe', {cpe_all});
end

function cg = coarse_grain(x, scale)
% 均值粗粒化：每个尺度窗口求平均
N = numel(x);
L = floor(N / scale);
if L < 1
    cg = [];
    return;
end
x = x(1:L*scale);
% 先 reshape 再对行求均值，可避免循环
cg = mean(reshape(x, scale, L), 1); % 1×L
cg = cg.'; % 转为列
end
