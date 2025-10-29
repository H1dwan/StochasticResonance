function dynamicsFunc = DynamicsFactory(systemType, systemParams)
% DynamicsFactory: 根据系统类型创建相应的动力学函数 dx/dt = -dU/dt + s + noise
% 输入参数:
%   systemType: 字符串，系统类型 ('CBSR', 'MUBSR')
%   systemParams: 结构体，包含系统参数
% 输出参数:
%   dynamicsFunc: 函数句柄，指向相应的动力学函数

switch upper(systemType)
    case 'CBSR'
        % 经典双稳态系统
        % 参数: a, b
        dynamicsFunc = @(x, ~) CBSR_Dynamics(x, ...
            systemParams.a, systemParams.b);
        
    case 'MUBSR'
        % 分段线性双稳态系统
        % 参数: a, b, c
        dynamicsFunc = @(x, ~) MUBSR_Dynamics(x, ...
            systemParams.a, systemParams.b, systemParams.c);
        
    case 'PLBSR'
        % PLBSR: 分段全线性双稳态系统
        % 参数: a, b, c
        dynamicsFunc = @(x, ~) PLBSR_Dynamics(x, ...
            systemParams.a, systemParams.b, systemParams.c);
        
    case 'PNBSR'
        % PNBSR: 分段非线性双稳态系统
        % 参数: a, b, c
        dynamicsFunc = @(x, ~) PNBSR_Dynamics(x, ...
            systemParams.a, systemParams.b, systemParams.c);
        
    case 'PEBSR'
        % PEBSR: 分段指数双稳态系统
        % 参数: a, b, c, l
        dynamicsFunc = @(x, ~) PEBSR_Dynamics(x, ...
            systemParams.a, systemParams.b, systemParams.c, systemParams.l);
        
    case 'MPBSR'
        % MPBSR: 分段线性多参数双稳态系统
        % 参数: p, u, c
        dynamicsFunc = @(x, ~) MPBSR_Dynamics(x, ...
            systemParams.p, systemParams.u, systemParams.c);
        
    otherwise
        % 用户自定义系统类型
        if isa(systemType, 'function_handle')
            dynamicsFunc = systemType;
        else
            error('不支持的系统类型: %s', systemType);
        end
end
end