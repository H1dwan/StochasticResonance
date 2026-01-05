function fig_handle = CreateThesisFigure(width_cm, height_cm, screen_scale)
% CreateThesisFigure  创建符合论文版式的 Figure（屏幕显示放大）
%
% 输入:
%   width_cm     - 图宽度（厘米），默认 8 cm
%   height_cm    - 图高度（厘米），默认 6 cm
%   screen_scale - 屏幕显示放大倍数，默认 2（只影响显示，不影响导出尺寸）
%
% 输出:
%   fig_handle   - 新建图句柄

if nargin < 1 || isempty(width_cm)
    width_cm = 8;
end
if nargin < 2 || isempty(height_cm)
    height_cm = 6;
end
if nargin < 3 || isempty(screen_scale)
    screen_scale = 2; % 放大显示，避免屏幕上过小
end

fig_handle = figure;

% 屏幕上的显示尺寸（仅显示放大，不改变导出比例）
set(fig_handle, 'Units', 'centimeters');
set(fig_handle, 'Position', [1, 1, width_cm * screen_scale, height_cm * screen_scale]);

% 背景色
set(fig_handle, 'Color', 'w');

% 打印/导出时的纸张设置（保证 PDF/PNG 比例一致）
set(fig_handle, 'PaperUnits', 'centimeters');
set(fig_handle, 'PaperPosition', [0, 0, width_cm, height_cm]);
set(fig_handle, 'PaperSize', [width_cm, height_cm]);
end
