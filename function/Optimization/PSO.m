function [fg, pg, Convergence_curve] = PSO(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% 粒子群算法（PSO）函数定义
% 输入参数：
% SearchAgents_no - 种群中搜索代理数量
% Max_iter - 最大迭代次数
% lb - 搜索空间下限
% ub - 搜索空间上限
% dim - 搜索空间维度
% fobj - 目标函数句柄
% 输出参数：
% fg - 全局最优解
% pg - 全局最优解位置
% Convergence_curve - 收敛曲线

% 初始化参数
Alpha_score = zeros(SearchAgents_no, 1); % 代理历史最优解的适应度
Alpha_pos = zeros(SearchAgents_no, dim); % 代理历史最优解的位置
Positions = zeros(SearchAgents_no, dim); % 代理位置
fitness = zeros(SearchAgents_no, 1); % 代理适应度
velocity = zeros(SearchAgents_no, dim); % 代理速度

pg = zeros(1, dim); % 全局最优解位置


vmax = abs(ub - lb); % 速度最大值
vmin = -vmax; % 速度最小值

% 初始化代理位置和速度
for i = 1:SearchAgents_no
    for j = 1:dim
        Positions(i,j) = lb(j) + rand() * (ub(j) - lb(j));
        velocity(i,j) = vmin(j) + rand() * (vmax(j) - vmin(j));  
    end
end

Convergence_curve = zeros(1, Max_iter); % 初始化收敛曲线

% 初始化代理得分和全局最优解位置
for i = 1:SearchAgents_no  
    fitness(i,1) = fobj(Positions(i,:));
    Alpha_score(i,1) = fitness(i,1);
    Alpha_pos(i,:) = Positions(i,:);
end

% 获取初始全局最优解位置和得分
[s, m] = min(Alpha_score(:,1));
fg = s;
pg(1,:) = Alpha_pos(m,:);

g = 0; % 循环计数器

% 开始迭代
while g < Max_iter
    w = 0.5 - (0.2 * g / Max_iter); % 更新惯性权重 0.5 0.2

    % 更新代理位置和速度
    for i = 1:SearchAgents_no
        for j = 1:dim     
            velocity(i,j) = w * velocity(i,j) + 2 * rand() * (Alpha_pos(i,j) - Positions(i,j)) + 2 * rand() * (pg(1,j) - Positions(i,j)); 
%             其中，velocity(i,j) 表示搜索点 $i$ 在第 $j$ 个维度上的速度，
%             Alpha_pos(i,j) 表示搜索点 $i$ 在第 $j$ 个维度上的个体最优位置，
%             pg(1,j) 表示全局最优位置在第 $j$ 个维度上的值，
%             Positions(i,j) 表示搜索点 $i$ 在第 $j$ 个维度上的位置，
%             w 是惯性权重。这段代码实现了上述速度更新公式中的所有项。
%             其中，2 * rand() * (Alpha_pos(i,j) - Positions(i,j)) 表示个体经验项，
%             2 * rand() * (pg(1,j) - Positions(i,j)) 表示社会经验项，
%             w * velocity(i,j) 表示惯性项。通过更新搜索点的速度，
%             可以让搜索点在下一步迭代中朝着更好的方向移动。
            if velocity(i,j) > vmax(j)
                velocity(i,j) = vmax(j) * rand();
            elseif velocity(i,j) < vmin(j)
                velocity(i,j) = vmin(j) * rand();
            end 

            Positions(i,j) = Positions(i,j) + velocity(i,j);

            % 处理越界情况
            if Positions(i,j) > ub(j)
                Positions(i,j) = ub(j) * rand();
            elseif Positions(i,j) < lb(j)
                Positions(i,j) = lb(j) * rand();
            end
        end
    end

    % 更新代理历史最优解的得分和位置
    for i = 1:SearchAgents_no
        fitness(i,1) = fobj(Positions(i,:));
        % 如果新得分优于旧得分，更新代理历史最优解的得分和位置
        if fitness(i,1) < Alpha_score(i,1) 
            Alpha_score(i,1) = fitness(i,1);
            Alpha_pos(i,:) = Positions(i,:);
        end
    end

    % 如果全局最优解得分有所改善，更新全局最优解位置
    if fg > min(Alpha_score(:,1))     
        [s, m] = min(Alpha_score(:,1));
        fg = s;
        pg(1,:) = Alpha_pos(m,:);
    end

    g = g + 1; % 更新迭代次数
    Convergence_curve(g) = fg; % 更新收敛曲线
end


