function ApplyThesisStyle(fig_handle, varargin)
% ApplyThesisStyle 将当前图表格式化为硕士论文标准风格 (v2.0)
%
% 输入:
%   fig_handle - 图形句柄 (通常用 gcf)
%   varargin   - 可选参数对，如 'FontSize', 10
%
% 更新日志:
%   v2.0: 字体改为Arial，刻度朝外，优化了Label的处理逻辑

% 解析参数
p = inputParser;
addRequired(p, 'fig_handle');
addParameter(p, 'FontSize', 10);   % 默认字号 10pt (适合单栏图)
addParameter(p, 'LineWidth', 1.0); % 坐标轴线宽 1.0pt
addParameter(p, 'GridOn', 'off');
parse(p, fig_handle, varargin{:});

font_name = 'Arial';  % 【更新】字体改为 Arial
font_size = p.Results.FontSize;
line_width = p.Results.LineWidth;
grid_state = p.Results.GridOn;

% 获取所有坐标轴
axes_handles = findall(fig_handle, 'type', 'axes');

for i = 1:length(axes_handles)
    ax = axes_handles(i);
    
    % --- 核心属性设置 ---
    set(ax, ...
        'FontName', font_name, ...
        'FontSize', font_size, ...
        'LineWidth', line_width, ...
        'Box', 'on', ...                   % 开启边框
        'XMinorTick', 'on', ...            % 开启次刻度
        'YMinorTick', 'on', ...
        'TickDir', 'out', ...              % 【更新】刻度线朝外
        'TickLength', [0.015 0.025]);      % 稍微缩短刻度线长度，配合朝外风格
    
    % --- 网格设置 ---
    if strcmpi(grid_state, 'on')
        grid(ax, 'on');
        set(ax, 'GridLineStyle', '--', 'GridAlpha', 0.3); % 降低网格透明度
    end
    
    % --- 标签强化 ---
    % 将 Label 和 Title 也统一为 Arial 并加粗
    set(ax.XLabel, 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    set(ax.YLabel, 'FontName', font_name, 'FontSize', font_size+1, 'FontWeight', 'bold');
    set(ax.Title,  'FontName', font_name, 'FontSize', font_size+2, 'FontWeight', 'bold');
    
    % --- 图例美化 ---
    if ~isempty(ax.Legend)
        set(ax.Legend, 'FontName', font_name, 'FontSize', font_size, 'Box', 'off');
    end
end

% 设置背景色为纯白
set(fig_handle, 'Color', 'w');

% 统一所有线条的粗细（防止默认 0.5 太细）
lines = findall(fig_handle, 'Type', 'line');
for k = 1:length(lines)
    if lines(k).LineWidth < 1.0
        lines(k).LineWidth = 1.5; % 默认数据曲线线宽
    end
end

end