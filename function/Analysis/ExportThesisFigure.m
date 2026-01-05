function ExportThesisFigure(fig_handle, file_base_name)
% ExportThesisFigure  按论文要求导出图像（PDF + SVG）
%
% 输入:
%   fig_handle      - 图句柄（可省略，默认 gcf）
%   file_base_name  - 文件基础名（不带扩展名）
%
% 说明:
%   将导出:
%       file_base_name.pdf  (矢量图)
%       file_base_name.svg  (矢量图)

if nargin < 1 || isempty(fig_handle)
    fig_handle = gcf;
end
if nargin < 2
    error('必须指定 file_base_name，例如 ''fig_SNR_D''');
end

% 导出为矢量 PDF
exportgraphics(fig_handle, [file_base_name, '.pdf'], ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

% 导出为高分辨率 svg
% 导出 SVG（R2021 用 print）
set(fig_handle,'Renderer','painters');
print(fig_handle, [file_base_name, '.svg'], '-dsvg');
end
