function SetThesisDefaultStyle()
% SetThesisDefaultStyle  设置硕士论文统一绘图风格
%
% 使用方法：
%   在任意脚本开头调用：
%       SetThesisDefaultStyle;

% ---------------- 全局字体与基本外观 ----------------
set(0, 'DefaultFigureColor',       'w');      % 图背景白色
set(0, 'DefaultAxesFontName',      'Arial');  % 坐标轴字体
set(0, 'DefaultTextFontName',      'Arial');  % 文本/标题字体
set(0, 'DefaultAxesFontSize',      10);       % 坐标轴刻度字号
set(0, 'DefaultAxesLabelFontSizeMultiplier', 1.2);
set(0, 'DefaultAxesTitleFontSizeMultiplier', 1.2);
set(0, 'DefaultTextFontSize',      12);       % 默认文本字号
set(0, 'DefaultAxesLineWidth',     1);        % 坐标轴线宽
set(0, 'DefaultLineLineWidth',     1.2);        % 曲线线宽
set(0, 'DefaultLineMarkerSize',    5);        % 标记大小

% 坐标轴习惯设置
set(0, 'DefaultAxesBox',           'on');     % 轴加框
set(0, 'DefaultAxesTickDir',       'out');    % 刻度朝外
set(0, 'DefaultAxesTickLabelInterpreter', 'none'); % 标签不使用 TeX
set(0, 'DefaultLegendBox',         'off');    % 默认图例无边框

% ---------------- 颜色顺序（统一配色） ----------------
color_order = [ ...
    % 0,      0,      0;        % 黑色
    % 0,      0.4470, 0.7410;   % 深蓝
    % 0.8500, 0.3250, 0.0980;   % 橙色
    % 0.4660, 0.6740, 0.1880;   % 绿色
    % 0.4940, 0.1840, 0.5560;   % 紫色
    0.9020, 0.2941, 0.2078; % 朱红
    0.2510, 0.3686, 0.6353; % 深蓝
    0.0000, 0.6196, 0.4510; % 墨绿
    0.5098, 0.2353, 0.5451; % 紫色
    0.9608, 0.5608, 0.2235; % 橙色
    0.2549, 0.7137, 0.7686; % 青色
    0.4941, 0.4941, 0.4941  % 灰色
    ];
set(0, 'DefaultAxesColorOrder', color_order);

% ---------------- 可选：使用 LaTeX 排版标签 ----------------
% 如果你经常用 LaTeX 公式，比如 '$SNR_{\mathrm{out}}$'
% 可以将下面三行取消注释：
set(0, 'DefaultTextInterpreter',        'latex');
set(0, 'DefaultAxesTickLabelInterpreter','latex');
set(0, 'DefaultLegendInterpreter',      'latex');
end
