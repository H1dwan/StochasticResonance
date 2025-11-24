function C = LZC(x)
% 计算信号的 Lempel-Ziv 复杂度 (LZC)
%
% Lempel-Ziv 复杂度是一种用于量化序列复杂性的指标，广泛应用于信号处理、
% 生物医学工程等领域。该算法通过计算序列中新子串出现的次数来衡量其复杂性。
%
% 输入参数:
%   x - 待分析的一维数值信号向量
%
% 输出参数:
%   C - Lempel-Ziv 复杂度值（未归一化）
%
% 参考文献:
%   [13] Lempel, A., & Ziv, J. (1976). On the complexity of finite sequences.
%   [14] Zhang, Y., et al. (2010). Detecting moving target in complex scenes.
%   [15] Wu, S., et al. (2013). Measuring complexity of short-term traffic.

% 1. 将信号二值化 (基于中位数) [13]
med = median(x);
s = x > med; % 转换为 0/1 序列

% 2. LZC 算法
% 初始化算法参数
n = length(s);
C = 1; % 复杂度计数器
i = 0; % 当前位置
k = 1; % 当前块长度
k_max = 1; % 最大块长度

% 主循环：查找新的子串模式
while i + k <= n
    % 检查 s(i+1...i+k) 是否在 s(1...i+k-1) 中出现过
    sub_str = s(i+1 : i+k);
    found = false;
    
    % 仅在 s(1...i+k-1) 中搜索
    for j = 1 : i
        if isequal(s(j : j+k-1), sub_str)
            found = true;
            break;
        end
    end
    
    if found
        % 如果找到了匹配的子串，则扩展当前块长度
        k = k + 1;
    else
        % 未找到匹配项，表示发现了一个新的子串模式
        % 增加复杂度计数，移动到下一个位置并重置块长度
        C = C + 1;
        i = i + k;
        k = 1;
    end
    
    % 更新最大块长度
    if k > k_max
        k_max = k;
    end
    
    % 边界处理
    if i + k > n
        if k > 1 % 如果 k > 1, 意味着最后一个块是重复的
            C = C + 1; % 最后一个新块
        end
        break;
    end
end

% 归一化 (可选, 但这里返回原始计数值)
C = C / (n / log2(n));
end