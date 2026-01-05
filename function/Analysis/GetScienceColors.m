function colors = GetScienceColors(palette_name)
% GetScienceColors 获取符合学术发表要求的配色方案
% 输入:
%   palette_name - 字符串，可选 'NPG' (Nature), 'AAAS' (Science), 'IEEE' (经典蓝红)
% 输出:
%   colors       - N×3 的 RGB 矩阵
%
% 示例:
%   cols = GetScienceColors('NPG');
%   plot(x, y, 'Color', cols(1,:));

if nargin < 1
    palette_name = 'NPG';
end

switch upper(palette_name)
    case 'NPG' 
        % Nature Publishing Group 风格 (推荐用于多条曲线对比)
        colors = [
            0.9020, 0.2941, 0.2078; % 朱红
            0.2510, 0.3686, 0.6353; % 深蓝
            0.0000, 0.6196, 0.4510; % 墨绿
            0.5098, 0.2353, 0.5451; % 紫色
            0.9608, 0.5608, 0.2235; % 橙色
            0.2549, 0.7137, 0.7686; % 青色
            0.4941, 0.4941, 0.4941  % 灰色
        ];
        
    case 'IEEE'
        % 经典的蓝红黑风格，适合双稳态势阱展示
        colors = [
            0.0000, 0.4470, 0.7410; % 标准蓝
            0.8500, 0.3250, 0.0980; % 标准红
            0.0000, 0.0000, 0.0000; % 纯黑
            0.4940, 0.1840, 0.5560  % 紫色
        ];
        
    otherwise
        % 默认 MATLAB 2019b+ 改进色系
        colors = lines(7);
end
end