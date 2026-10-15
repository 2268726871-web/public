function experiment1_equilibrium_validation_restart()
%EXPERIMENT1_EQUILIBRIUM_VALIDATION_RESTART
% 论文仿真实验E1：六仓24 h非线性平衡点验证。
%
% 验证条件：
%   1. 室外环境和生物源保持标称值；
%   2. 六个仓从各自阶段目标状态开始；
%   3. 控制输入固定为共享参数文件计算的全精度平衡输入；
%   4. 不启用反馈、协同、事件触发和安全强制触发。

clc;
close all;

%% 1. 定位统一参数文件
scriptDir = fileparts(mfilename('fullpath'));
defaultCodeDir = '/Users/bystander/Desktop/蘑菇仓/代码库';

if exist(fullfile(scriptDir,'mushroom_system_parameters.m'),'file') == 2
    codeDir = scriptDir;
elseif exist(fullfile(defaultCodeDir,'mushroom_system_parameters.m'),'file') == 2
    codeDir = defaultCodeDir;
else
    error('未找到mushroom_system_parameters.m，请检查代码库位置。');
end

addpath(codeDir);
cfg = mushroom_system_parameters();

%% 2. 只检查本实验需要的参数
assert(cfg.room.nRoom == 6,'本实验要求六个蘑菇仓。');
assert(cfg.stage.nStage == 4,'本实验要求G1-G4四个阶段。');
assert(all(isfinite(cfg.equilibrium.u0(:))), ...
    '全精度平衡输入中存在NaN或Inf。');
assert(all(isfinite(cfg.audit.balanceResidual(:))), ...
    '平衡方程残差中存在NaN或Inf。');

dt_s = cfg.time.dt_s;
tEnd_s = cfg.time.equilibriumDuration_s;
t_s = (0:dt_s:tEnd_s)';
nStep = numel(t_s);
nRoom = cfg.room.nRoom;
roomStage = cfg.room.stageIndex;

% 全精度平衡验证容差。
tolT_C = 1e-6;
tolW_kg_kgda = 1e-10;
tolCO2_ppm = 1e-3;

%% 3. 标称室外状态
out.T_C = cfg.outdoor.nominal.T_C;
out.w_kg_kgda = cfg.outdoor.nominal.w_kg_kgda;
out.cCO2_kg_m3 = cfg.outdoor.nominal.cCO2_kg_m3;

%% 4. 六仓从各阶段目标状态开始
T_C = zeros(nStep,nRoom);
w_kg_kgda = zeros(nStep,nRoom);
cCO2_kg_m3 = zeros(nStep,nRoom);

for i = 1:nRoom
    g = roomStage(i);
    T_C(1,i) = cfg.stage.TTarget_C(g);
    w_kg_kgda(1,i) = cfg.stage.wTarget_kg_kgda(g);
    cCO2_kg_m3(1,i) = cfg.stage.cCO2Target_kg_m3(g);
end

%% 5. 初始平衡残差
initialResidual = zeros(nRoom,3);

for i = 1:nRoom
    g = roomStage(i);
    x0 = [T_C(1,i);w_kg_kgda(1,i);cCO2_kg_m3(1,i)];
    u0 = cfg.equilibrium.u0(g,:)';
    initialResidual(i,:) = room_rhs(x0,u0,g,out,cfg)';
end

%% 6. 使用RK4积分24 h
for k = 1:(nStep-1)
    for i = 1:nRoom
        g = roomStage(i);
        xk = [T_C(k,i);w_kg_kgda(k,i);cCO2_kg_m3(k,i)];
        u0 = cfg.equilibrium.u0(g,:)';

        k1 = room_rhs(xk,u0,g,out,cfg);
        k2 = room_rhs(xk+0.5*dt_s*k1,u0,g,out,cfg);
        k3 = room_rhs(xk+0.5*dt_s*k2,u0,g,out,cfg);
        k4 = room_rhs(xk+dt_s*k3,u0,g,out,cfg);

        xNext = xk + dt_s*(k1+2*k2+2*k3+k4)/6;
        T_C(k+1,i) = xNext(1);
        w_kg_kgda(k+1,i) = xNext(2);
        cCO2_kg_m3(k+1,i) = xNext(3);
    end
end

%% 7. 转换CO2并计算漂移
CO2_ppm = zeros(nStep,nRoom);

for i = 1:nRoom
    for k = 1:nStep
        CO2_ppm(k,i) = mass_to_ppm_local( ...
            cCO2_kg_m3(k,i),T_C(k,i), ...
            cfg.physical.PIndoorRef_Pa,cfg.physical);
    end
end

targetT_C = cfg.stage.TTarget_C(roomStage);
targetW_kg_kgda = cfg.stage.wTarget_kg_kgda(roomStage);
targetCO2_ppm = cfg.stage.CO2Target_ppm(roomStage);

dT_C = bsxfun(@minus,T_C,targetT_C');
dW_kg_kgda = bsxfun(@minus,w_kg_kgda,targetW_kg_kgda');
dCO2_ppm = bsxfun(@minus,CO2_ppm,targetCO2_ppm');

maxAbsT_C = max(abs(dT_C),[],1)';
maxAbsW_kg_kgda = max(abs(dW_kg_kgda),[],1)';
maxAbsCO2_ppm = max(abs(dCO2_ppm),[],1)';

finalDT_C = dT_C(end,:)';
finalDW_kg_kgda = dW_kg_kgda(end,:)';
finalDCO2_ppm = dCO2_ppm(end,:)';

passFlag = maxAbsT_C <= tolT_C ...
    & maxAbsW_kg_kgda <= tolW_kg_kgda ...
    & maxAbsCO2_ppm <= tolCO2_ppm;

%% 8. 汇总表
stageTable = table( ...
    cfg.stage.names, ...
    cfg.equilibrium.Q0_W, ...
    cfg.equilibrium.mWater0_kg_s, ...
    3600*cfg.equilibrium.mWater0_kg_s, ...
    cfg.equilibrium.F0_m3_s, ...
    3600*cfg.equilibrium.F0_m3_s, ...
    cfg.audit.balanceResidual(:,1), ...
    cfg.audit.balanceResidual(:,2), ...
    cfg.audit.balanceResidual(:,3), ...
    'VariableNames',{ ...
    'Stage','Q0_W','mWater0_kg_s','mWater0_kg_h', ...
    'F0_m3_s','F0_m3_h','ResidualT','ResidualW','ResidualCO2'});

roomTable = table( ...
    cfg.room.names,cfg.room.stageName, ...
    targetT_C,targetW_kg_kgda,targetCO2_ppm, ...
    finalDT_C,finalDW_kg_kgda,finalDCO2_ppm, ...
    maxAbsT_C,maxAbsW_kg_kgda,maxAbsCO2_ppm,passFlag, ...
    'VariableNames',{ ...
    'Room','Stage','TargetT_C','TargetW_kg_kgda','TargetCO2_ppm', ...
    'FinalDT_C','FinalDW_kg_kgda','FinalDCO2_ppm', ...
    'MaxAbsDT_C','MaxAbsDW_kg_kgda','MaxAbsDCO2_ppm','Pass'});

residualTable = table( ...
    cfg.room.names,cfg.room.stageName, ...
    initialResidual(:,1),initialResidual(:,2),initialResidual(:,3), ...
    'VariableNames',{ ...
    'Room','Stage','dT_dt_C_s','dw_dt_kg_kgda_s','dc_dt_kg_m3_s'});

%% 9. 保存到全新的论文重跑目录
resultDir = fullfile(codeDir,'results','paper_rerun','E1_equilibrium');

if exist(resultDir,'dir') ~= 7
    mkdir(resultDir);
end

writetable(stageTable,fullfile(resultDir,'stage_equilibrium_inputs.csv'));
writetable(roomTable,fullfile(resultDir,'room_equilibrium_validation.csv'));
writetable(residualTable,fullfile(resultDir,'initial_balance_residual.csv'));

save(fullfile(resultDir,'E1_equilibrium_validation.mat'), ...
    'cfg','t_s','T_C','w_kg_kgda','cCO2_kg_m3','CO2_ppm', ...
    'dT_C','dW_kg_kgda','dCO2_ppm','stageTable', ...
    'roomTable','residualTable');

%% 10. 绘制六仓平衡漂移图
time_h = t_s/3600;
lineColor = lines(nRoom);
fig = figure('Color','w','Position',[100 100 1180 800]);

subplot(3,1,1);
hold on;
for i = 1:nRoom
    plot(time_h,dT_C(:,i),'LineWidth',1.4,'Color',lineColor(i,:));
end
grid on;
ylabel('温度漂移/°C');
title('E1 六仓24 h非线性平衡点验证');
legend(cfg.room.names,'Location','eastoutside');

subplot(3,1,2);
hold on;
for i = 1:nRoom
    plot(time_h,1e6*dW_kg_kgda(:,i), ...
        'LineWidth',1.4,'Color',lineColor(i,:));
end
grid on;
ylabel('含湿量漂移/(10^{-6} kg/kg)');

subplot(3,1,3);
hold on;
for i = 1:nRoom
    plot(time_h,dCO2_ppm(:,i),'LineWidth',1.4,'Color',lineColor(i,:));
end
grid on;
xlabel('时间/h');
ylabel('CO₂漂移/ppm');

saveas(fig,fullfile(resultDir,'E1_equilibrium_drift.png'));

%% 11. 命令行结论
fprintf('\n============================================================\n');
fprintf('E1 六仓24 h非线性平衡点验证完成\n');
fprintf('参数来源：mushroom_system_parameters.m\n');
fprintf('平衡输入：程序计算的全精度值\n');
fprintf('结果目录：%s\n',resultDir);
fprintf('最大温度漂移：%.15g °C\n',max(maxAbsT_C));
fprintf('最大含湿量漂移：%.15g kg/kg\n',max(maxAbsW_kg_kgda));
fprintf('最大CO₂漂移：%.15g ppm\n',max(maxAbsCO2_ppm));
fprintf('============================================================\n');

disp(stageTable);
disp(roomTable);

if all(passFlag)
    fprintf('结论：PASS，六个仓均通过24 h平衡点验证。\n\n');
else
    error('结论：FAIL，至少一个仓超过平衡点验证容差。');
end
end

function dx = room_rhs(x,u,g,out,cfg)
T_C = x(1);
w_kg_kgda = x(2);
cCO2_kg_m3 = x(3);

Q_W = u(1);
mWater_kg_s = u(2);
F_m3_s = u(3);

dT_C_s = ( ...
    cfg.physical.UA_W_K*(out.T_C-T_C) ...
    + cfg.physical.rhoAir_kg_m3*cfg.physical.cpAir_J_kgK ...
      *F_m3_s*(out.T_C-T_C) ...
    + Q_W + cfg.bio.Q_W(g)) / cfg.physical.CT_J_K;

dw_kg_kgda_s = ( ...
    cfg.physical.rhoAir_kg_m3*F_m3_s*(out.w_kg_kgda-w_kg_kgda) ...
    + mWater_kg_s + cfg.bio.mWater_kg_s(g)) / cfg.stage.Mda_kg(g);

dcCO2_kg_m3_s = ( ...
    F_m3_s*(out.cCO2_kg_m3-cCO2_kg_m3) ...
    + cfg.bio.mCO2_kg_s(g)) / cfg.physical.V_m3;

dx = [dT_C_s;dw_kg_kgda_s;dcCO2_kg_m3_s];
end

function ppm = mass_to_ppm_local(c,T_C,P_Pa,physical)
ppm = c*physical.R_J_molK*(T_C+273.15) ...
    /(P_Pa*physical.MCO2_kg_mol)*1e6;
end
