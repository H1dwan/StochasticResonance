function [best_fit, best_theta, history] = PGPSO(theta_min, theta_max, opts, evaluator)
% PGPSO 基于物理知识引导的自适应粒子群优化算法搜索最优参数
%   [best_fit, best_theta, history] = PGPSO(theta_min, theta_max, opts, evaluator)
%   在给定参数上下界与评估器 evaluator 下执行盲指标引导的粒子群优化（PSO），自适应调整惯性权重
%   与加速度系数以寻找最佳参数组合。
%
% 输入参数
%   theta_min   - 参数向量下界（列向量或行向量）。
%   theta_max   - 参数向量上界（列向量或行向量）。
%   opts        - 结构体，包含以下字段：
%       .rand_seed     - 随机种子（可选）。
%       .num_particles - 粒子数量。
%       .max_iter      - 最大迭代次数。
%       .w_min, .w_max - 惯性权重范围。
%       .c1_min, .c1_max - 认知加速度系数范围。
%       .c2_min, .c2_max - 社会加速度系数范围
%       .display       - 是否显示迭代信息（布尔值）。
%       .goal          - 优化方向，'max'（默认）或 'min'。
%   evaluator   - 函数句柄，接受参数向量，返回目标适应度值、共振度和调试信息（可选）。
%
% 输出参数
%   best_fit    - 最佳适应度值。
%   best_theta  - 最佳参数向量。
%   history     - 结构体，包含以下字段：
%       .gbest_fit - 每次迭代的全局最佳适应度值。
%       .mean_fit  - 每次迭代的平均适应度值。

%% 1. 随机种子设置 ---------------------------------------------------
if isfield(opts, 'rand_seed')
    rng(opts.rand_seed);
end
rng(1);  % 固定随机种子以便复现
%% 2. 读取 PSO 与盲指标相关参数 ------------------------------------
num_particles = opts.num_particles;
max_iter      = opts.max_iter;

w_min   = opts.w_min;
w_max   = opts.w_max;
c1_min  = opts.c1_min;
c1_max  = opts.c1_max;
c2_min  = opts.c2_min;
c2_max  = opts.c2_max;

do_display   = opts.display;
is_minimize  = isfield(opts, 'goal') && strcmpi(opts.goal, 'min');
if is_minimize
    better_than = @(new, old) new < old;
    select_best = @min;
else
    better_than = @(new, old) new > old;
    select_best = @max;
end

%% 3. 粒子初始化 -----------------------------------------------------
theta_min = theta_min(:)';  % 确保为行向量
theta_max = theta_max(:)';
dim_theta = numel(theta_min);

% 位置初始化：在参数区间内均匀随机
theta_particles = rand(num_particles, dim_theta) ...
    .* (theta_max - theta_min) + theta_min;

% 初始速度设为 0
vel_particles   = zeros(num_particles, dim_theta);

% 个体最优位置与适应度
pbest_pos = theta_particles;
if is_minimize
    pbest_fit = inf(num_particles, 1);
    gbest_fit = inf;
else
    pbest_fit = -inf(num_particles, 1);
    gbest_fit = -inf;
end

% 全局最优
gbest_pos = theta_particles(1,:);

% 历史记录
history.gbest_fit = zeros(max_iter, 1);
history.mean_fit  = zeros(max_iter, 1);
history.mean_eta    = zeros(max_iter, 1);   % 新增：群体平均 eta
history.best_eta    = zeros(max_iter, 1);   % 新增：当前全局最优粒子的 eta
history.corr_etaFit = zeros(max_iter, 1);   % 新增：eta 与 J 相关系数

%% 4. 初始评估 -------------------------------------------------------
fit_vals = zeros(num_particles, 1);
eta_vals = zeros(num_particles, 1);

for idx = 1:num_particles
    [fit_vals(idx), eta_vals(idx)] = evaluator(theta_particles(idx,:));
    
    pbest_fit(idx)   = fit_vals(idx);
    pbest_pos(idx,:) = theta_particles(idx,:);
    
    if better_than(fit_vals(idx), gbest_fit)
        gbest_fit = fit_vals(idx);
        gbest_pos = theta_particles(idx,:);
    end
end

history.gbest_fit(1) = gbest_fit;
history.mean_fit(1)  = mean(fit_vals);
history.mean_eta(1)  = mean(eta_vals);

if num_particles >= 2
    history.corr_etaFit(1) = corr(eta_vals, fit_vals);
else
    history.corr_etaFit(1) = NaN;
end

% 全局最优粒子的 eta
[~, idx_best_init]   = select_best(fit_vals);
history.best_eta(1)  = eta_vals(idx_best_init);

if do_display
    fprintf('Iter %3d, best J = %.4f, mean J = %.4f\n', ...
        1, history.gbest_fit(1), history.mean_fit(1));
end

%% 5. PSO 主循环 -----------------------------------------------------
for iter = 2:max_iter
    
    for idx = 1:num_particles
        
        % -------- 由盲共振度 eta 自适应 w, c1, c2 ---------------
        eta_i = eta_vals(idx);             % (0,1)，越大认为越接近共振
        
        w_i  = w_min + (w_max - w_min) * (1 - eta_i);
        c1_i = c1_min + (c1_max - c1_min) * eta_i;
        c2_i = c2_min + (c2_max - c2_min) * eta_i;
        
        r1 = rand;
        r2 = rand;
        
        % -------- 速度更新 -----------------------------------------
        vel_particles(idx,:) = ...
            w_i  * vel_particles(idx,:) ...
            + c1_i * r1 * (pbest_pos(idx,:) - theta_particles(idx,:)) ...
            + c2_i * r2 * (gbest_pos        - theta_particles(idx,:));
        
        % -------- 位置更新 + 边界裁剪 -------------------------------
        theta_particles(idx,:) = theta_particles(idx,:) + vel_particles(idx,:);
        
        theta_particles(idx,:) = max(theta_particles(idx,:), theta_min);
        theta_particles(idx,:) = min(theta_particles(idx,:), theta_max);
        
        % -------- 重新评估盲适应度 ---------------------------------
        [fit_vals(idx), eta_vals(idx)] = evaluator(theta_particles(idx,:));
        
        % -------- 更新个体最优 -------------------------------------
        if better_than(fit_vals(idx), pbest_fit(idx))
            pbest_fit(idx)   = fit_vals(idx);
            pbest_pos(idx,:) = theta_particles(idx,:);
        end
    end
    
    % -------- 更新全局最优 -----------------------------------------
    [best_now, idx_best_now] = select_best(fit_vals);
    if better_than(best_now, gbest_fit)
        gbest_fit = best_now;
        gbest_pos = theta_particles(idx_best_now,:);
    end
    
    history.gbest_fit(iter) = gbest_fit;
    history.mean_fit(iter)  = mean(fit_vals);
    
    history.mean_eta(iter)  = mean(eta_vals);
    
    if num_particles >= 2
        history.corr_etaFit(iter) = corr(eta_vals, fit_vals);
    else
        history.corr_etaFit(iter) = NaN;
    end
    
    [~, idx_best_now]    = select_best(fit_vals);
    history.best_eta(iter) = eta_vals(idx_best_now);
    
    if do_display
        fprintf('Iter %3d, best J = %.4f, mean J = %.4f, mean eta = %.3f, corr(eta,J) = %.3f\n', ...
            iter, history.gbest_fit(iter), history.mean_fit(iter), ...
            history.mean_eta(iter), history.corr_etaFit(iter));
    end
end

best_theta = gbest_pos;
best_fit   = gbest_fit;

end