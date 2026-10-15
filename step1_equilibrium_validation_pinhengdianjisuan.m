function step1_equilibrium_validation
%STEP1_EQUILIBRIUM_VALIDATION 六个蘑菇仓的24 h平衡点验证。
%
% 本实验只验证物理模型和平衡输入是否一致：
%   1) 室外状态和生物扰动保持为标称值；
%   2) 每个仓从所属阶段的目标状态开始；
%   3) 控制输入固定为对应阶段的平衡输入；
%   4) 不启用反馈控制、仓间协同和事件触发。
%
% 运行方法：
%   1) 在 MATLAB 中把当前目录切换到本文件所在目录；
%   2) 输入：step1_equilibrium_validation
%
% 输出位于：simulation/results/step1_equilibrium/

clc;
close all;

%% 1. 仿真设置
dt = 60;                    % 数值积分步长，s
tEnd = 24 * 3600;           % 平衡点验证时长，s
t = (0:dt:tEnd)';
nStep = numel(t);
nRoom = 6;

% 这些只是数值验证容差，不是文献给出的生长安全范围。
tolT = 0.01;                % 最大允许温度漂移，degC
tolW = 1.0e-6;              % 最大允许含湿量漂移，kg/kg干空气
tolCppm = 1.0;              % 最大允许CO2漂移，ppm

% true：使用论文表格中保留位数后的平衡输入，检验四舍五入影响。
% false：使用程序重新计算的全精度平衡输入，理论漂移应接近机器精度。
useRoundedTableInputs = true;

%% 2. 物理常数和仓体参数
par.rhoAir = 1.2;           % 湿空气密度，kg/m^3
par.cpAir = 1006;           % 空气定压比热，J/(kg K)
par.UA = 31.5;              % 围护结构综合传热系数乘面积，W/K
par.CT = 6.85e6;            % 有效热容量，J/K
par.V = 104;                % 仓内容积，m^3

% CO2的ppm与质量浓度换算常数。
par.R = 8.314462618;        % 通用气体常数，J/(mol K)
par.MCO2 = 44.01e-3;        % CO2摩尔质量，kg/mol
par.PIndoorRef = 101325;    % 室内ppm换算参考压力，Pa

%% 3. 标称室外状态
out.T = 22.34160274;        % 室外温度，degC
out.w = 0.0146041409999;    % 室外含湿量，kg/kg干空气
out.P = 99249.5605;         % 南宁2023逐时气压平均值，Pa
out.Cppm = 419.35;          % 室外CO2背景值，ppm
out.c = ppm_to_mass(out.Cppm, out.T, out.P, par);

%% 4. G1-G4阶段目标和生物扰动
stageNames = {'G1'; 'G2'; 'G3'; 'G4'};
T0 = [23; 24; 16; 16];
w0 = [0.0168318; 0.0179081; 0.0096345; 0.0096345];
Cppm0 = [12500; 12500; 900; 900];
Mda = [120.7; 120.1; 125.0; 125.0];

Qbio = [170; 246; 246; 189];                       % W
mbioW = [7.87e-5; 8.36e-5; 1.53e-4; 1.53e-4];     % kg/s
mbioC = [4.73e-5; 6.83e-5; 6.83e-5; 5.25e-5];     % kg/s

% 室内目标ppm按既定室内参考压力换算为CO2质量浓度。
c0 = zeros(4, 1);
for g = 1:4
    c0(g) = ppm_to_mass(Cppm0(g), T0(g), par.PIndoorRef, par);
end

%% 5. 由三个平衡方程重新计算全精度平衡输入
% CO2平衡：0 = F0*(cout-c0) + mbioC
F0Exact = mbioC ./ (c0 - out.c);                    % m^3/s

% 含湿量平衡：0 = rho*F0*(wout-w0) + mw0 + mbioW
mw0Exact = par.rhoAir .* F0Exact .* (w0 - out.w) - mbioW; % kg/s

% 温度平衡：0 = UA*(Tout-T0) + rho*cp*F0*(Tout-T0)
%                  + Q0 + Qbio
Q0Exact = (par.UA + par.rhoAir .* par.cpAir .* F0Exact) ...
          .* (T0 - out.T) - Qbio;                   % W

% 已确定的表格值（经过显示位数四舍五入）。
Q0Rounded = [-147.543; -187.493; -1011.754; -823.821];
mw0Rounded = [-7.29243e-5; -7.11874e-5; -5.93899e-4; -4.91905e-4];
F0Rounded = [0.002160593; 0.003130736; 0.073932066; 0.056829187];

if useRoundedTableInputs
    Q0Used = Q0Rounded;
    mw0Used = mw0Rounded;
    F0Used = F0Rounded;
    inputMode = '表格四舍五入值';
else
    Q0Used = Q0Exact;
    mw0Used = mw0Exact;
    F0Used = F0Exact;
    inputMode = '程序计算全精度值';
end

%% 6. 六仓阶段配置
% 1、2号为G1，3号为G2，4、5号为G3，6号为G4。
roomNames = {'1号仓'; '2号仓'; '3号仓'; '4号仓'; '5号仓'; '6号仓'};
roomStage = [1; 1; 2; 3; 3; 4];

% 注意：平衡点验证从阶段目标开始，而不是从后续控制实验的偏差初值开始。
T = zeros(nStep, nRoom);
w = zeros(nStep, nRoom);
c = zeros(nStep, nRoom);

for i = 1:nRoom
    g = roomStage(i);
    T(1, i) = T0(g);
    w(1, i) = w0(g);
    c(1, i) = c0(g);
end

%% 7. 使用四阶Runge-Kutta法积分24 h
initialResidual = zeros(nRoom, 3);

for i = 1:nRoom
    g = roomStage(i);
    x0Room = [T(1, i); w(1, i); c(1, i)];
    uRoom = [Q0Used(g); mw0Used(g); F0Used(g)];
    bioRoom = [Qbio(g); mbioW(g); mbioC(g)];
    initialResidual(i, :) = room_rhs(x0Room, uRoom, bioRoom, Mda(g), out, par)';
end

for k = 1:(nStep - 1)
    for i = 1:nRoom
        g = roomStage(i);
        xk = [T(k, i); w(k, i); c(k, i)];
        uRoom = [Q0Used(g); mw0Used(g); F0Used(g)];
        bioRoom = [Qbio(g); mbioW(g); mbioC(g)];

        k1 = room_rhs(xk, uRoom, bioRoom, Mda(g), out, par);
        k2 = room_rhs(xk + 0.5 * dt * k1, uRoom, bioRoom, Mda(g), out, par);
        k3 = room_rhs(xk + 0.5 * dt * k2, uRoom, bioRoom, Mda(g), out, par);
        k4 = room_rhs(xk + dt * k3, uRoom, bioRoom, Mda(g), out, par);
        xNext = xk + dt / 6 * (k1 + 2*k2 + 2*k3 + k4);

        T(k + 1, i) = xNext(1);
        w(k + 1, i) = xNext(2);
        c(k + 1, i) = xNext(3);
    end
end

%% 8. 将CO2质量浓度转换回ppm，并计算漂移量
Cppm = zeros(nStep, nRoom);
for i = 1:nRoom
    for k = 1:nStep
        Cppm(k, i) = mass_to_ppm(c(k, i), T(k, i), par.PIndoorRef, par);
    end
end

targetT = zeros(1, nRoom);
targetW = zeros(1, nRoom);
targetCppm = zeros(1, nRoom);
for i = 1:nRoom
    g = roomStage(i);
    targetT(i) = T0(g);
    targetW(i) = w0(g);
    targetCppm(i) = Cppm0(g);
end

dT = T - targetT;
dW = w - targetW;
dCppm = Cppm - targetCppm;

maxAbsT = max(abs(dT), [], 1)';
maxAbsW = max(abs(dW), [], 1)';
maxAbsCppm = max(abs(dCppm), [], 1)';
finalDT = dT(end, :)';
finalDW = dW(end, :)';
finalDCppm = dCppm(end, :)';
passFlag = maxAbsT <= tolT & maxAbsW <= tolW & maxAbsCppm <= tolCppm;

%% 9. 保存结果
thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(thisDir, 'results', 'step1_equilibrium');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

stageInputTable = table(stageNames, Q0Exact, Q0Rounded, ...
    mw0Exact, mw0Rounded, mw0Exact*3600, F0Exact, F0Rounded, F0Exact*3600, ...
    'VariableNames', {'Stage', 'Q0_exact_W', 'Q0_rounded_W', ...
    'mw0_exact_kg_s', 'mw0_rounded_kg_s', 'mw0_exact_kg_h', ...
    'F0_exact_m3_s', 'F0_rounded_m3_s', 'F0_exact_m3_h'});
writetable(stageInputTable, fullfile(resultDir, 'equilibrium_inputs.csv'));

roomStageNames = stageNames(roomStage);
summaryTable = table(roomNames, roomStageNames, targetT', targetW', targetCppm', ...
    finalDT, finalDW, finalDCppm, maxAbsT, maxAbsW, maxAbsCppm, passFlag, ...
    'VariableNames', {'Room', 'Stage', 'Target_T_C', 'Target_w_kg_kgda', ...
    'Target_CO2_ppm', 'Final_dT_C', 'Final_dw_kg_kgda', 'Final_dCO2_ppm', ...
    'MaxAbs_dT_C', 'MaxAbs_dw_kg_kgda', 'MaxAbs_dCO2_ppm', 'Pass'});
writetable(summaryTable, fullfile(resultDir, 'equilibrium_validation_summary.csv'));

residualTable = table(roomNames, roomStageNames, initialResidual(:,1), ...
    initialResidual(:,2), initialResidual(:,3), ...
    'VariableNames', {'Room', 'Stage', 'dT_dt_C_s', ...
    'dw_dt_kg_kgda_s', 'dc_dt_kg_m3_s'});
writetable(residualTable, fullfile(resultDir, 'initial_balance_residual.csv'));

save(fullfile(resultDir, 'equilibrium_validation.mat'), ...
    't', 'T', 'w', 'c', 'Cppm', 'dT', 'dW', 'dCppm', ...
    'Q0Exact', 'mw0Exact', 'F0Exact', 'summaryTable');

%% 10. 绘制平衡点漂移曲线
timeHour = t / 3600;
lineColors = lines(nRoom);

fig = figure('Color', 'w', 'Position', [100, 100, 1100, 780]);

subplot(3,1,1);
hold on;
for i = 1:nRoom
    plot(timeHour, dT(:,i), 'LineWidth', 1.4, 'Color', lineColors(i,:));
end
grid on;
ylabel('\DeltaT (^{\circ}C)');
title(['24 h平衡点验证：', inputMode]);
legend(roomNames, 'Location', 'eastoutside');

subplot(3,1,2);
hold on;
for i = 1:nRoom
    plot(timeHour, dW(:,i)*1e6, 'LineWidth', 1.4, 'Color', lineColors(i,:));
end
grid on;
ylabel('\Deltaw (10^{-6} kg/kg_{da})');

subplot(3,1,3);
hold on;
for i = 1:nRoom
    plot(timeHour, dCppm(:,i), 'LineWidth', 1.4, 'Color', lineColors(i,:));
end
grid on;
xlabel('时间 (h)');
ylabel('\DeltaCO_2 (ppm)');

saveas(fig, fullfile(resultDir, 'equilibrium_drift.png'));

%% 11. 命令行汇总
fprintf('\n============================================================\n');
fprintf('六仓24 h平衡点验证完成\n');
fprintf('平衡输入模式：%s\n', inputMode);
fprintf('结果目录：%s\n', resultDir);
fprintf('验证容差：|dT|<=%.3g degC, |dw|<=%.3g kg/kgda, |dCO2|<=%.3g ppm\n', ...
    tolT, tolW, tolCppm);
fprintf('============================================================\n');
disp(stageInputTable);
disp(summaryTable);

if all(passFlag)
    fprintf('结论：PASS，六个仓均通过平衡点验证。\n\n');
else
    warning('结论：FAIL，至少一个仓超过数值验证容差，请检查结果表。');
end
end

function dx = room_rhs(x, u, bio, Mda, out, par)
% 单仓温度-含湿量-CO2连续非线性模型。
T = x(1);
w = x(2);
c = x(3);

QT = u(1);
mw = u(2);
Fv = u(3);

Qbio = bio(1);
mbioW = bio(2);
mbioC = bio(3);

dT = (par.UA*(out.T - T) ...
      + par.rhoAir*par.cpAir*Fv*(out.T - T) ...
      + QT + Qbio) / par.CT;

dw = (par.rhoAir*Fv*(out.w - w) + mw + mbioW) / Mda;

dc = (Fv*(out.c - c) + mbioC) / par.V;

dx = [dT; dw; dc];
end

function c = ppm_to_mass(Cppm, TdegC, pressurePa, par)
% 理想气体关系：ppm -> kg/m^3。
c = Cppm * 1e-6 * pressurePa * par.MCO2 ...
    / (par.R * (TdegC + 273.15));
end

function Cppm = mass_to_ppm(c, TdegC, pressurePa, par)
% 理想气体关系：kg/m^3 -> ppm。
Cppm = c * par.R * (TdegC + 273.15) ...
    / (pressurePa * par.MCO2) * 1e6;
end
