function experiment2_continuous_linearization_validation()
%EXPERIMENT2_CONTINUOUS_LINEARIZATION_VALIDATION
% 论文仿真实验E2：G1-G4四阶段连续局部线性化核验。
%
% 本实验验证：
%   1. 解析连续矩阵A、B、Cout均为有限数；
%   2. 解析矩阵与非线性模型的中心差分Jacobian一致；
%   3. 平衡点附近的一阶线性化余项随扰动尺度按二阶规律缩小；
%   4. 四阶段连续A矩阵的全部特征值实部均小于0；
%   5. 物理状态到归一化输出的变换矩阵M均可逆且条件数合理。
%
% 说明：本实验只验证连续局部模型，不进行300 s零阶保持离散化。

clc;
close all;
format long g;

%% 1. 定位统一参数文件和连续矩阵生成函数
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

assert(exist('build_continuous_matrices','file') == 2, ...
    '未找到build_continuous_matrices.m。');

cfg = mushroom_system_parameters();
cont = build_continuous_matrices(cfg);

assert(cfg.stage.nStage == 4,'本实验要求G1-G4四个阶段。');
assert(size(cont.A,1) == 3 && size(cont.A,2) == 3, ...
    '连续状态矩阵A的尺寸应为3×3×4。');
assert(size(cont.B,1) == 3 && size(cont.B,2) == 3, ...
    '连续输入矩阵B的尺寸应为3×3×4。');

%% 2. 数值核验设置
nStage = cfg.stage.nStage;
nState = 3;
nInput = 3;

% 中心差分步长，依次对应T、w、c和Q、mWater、F。
stepX = [1e-2;1e-5;1e-5];
stepU = [1;1e-6;1e-5];

% 允许的解析矩阵与数值Jacobian相对误差。
tolRelA = 1e-6;
tolRelB = 1e-6;
tolRelC = 1e-6;
tolConditionM = 1e8;
tolOrder = 1.8;

% 局部线性化扰动尺度；一阶线性化余项理论上与alpha^2同阶。
alpha = [1;0.5;0.25;0.125];
nAlpha = numel(alpha);

% 在归一化状态空间中设置统一扰动方向，保证四阶段比较口径一致。
zDirection = [0.20;-0.15;0.10];

% 输入扰动基准：[热功率W；加/除湿kg/s；通风m3/s]。
uDirection = [100;1e-5;1e-3];

out.T_C = cfg.outdoor.nominal.T_C;
out.w_kg_kgda = cfg.outdoor.nominal.w_kg_kgda;
out.cCO2_kg_m3 = cfg.outdoor.nominal.cCO2_kg_m3;

%% 3. 四阶段解析矩阵与数值Jacobian对照
relErrorA = zeros(nStage,1);
relErrorB = zeros(nStage,1);
relErrorC = zeros(nStage,1);
maxRealEigA = zeros(nStage,1);
spectralAbscissa = zeros(nStage,1);
conditionM = zeros(nStage,1);
localRemainder = zeros(nStage,nAlpha);
observedOrder = zeros(nStage,nAlpha-1);

numericA = zeros(nState,nState,nStage);
numericB = zeros(nState,nInput,nStage);
numericCout = zeros(nState,nState,nStage);

for g = 1:nStage
    x0 = [ ...
        cfg.stage.TTarget_C(g); ...
        cfg.stage.wTarget_kg_kgda(g); ...
        cfg.stage.cCO2Target_kg_m3(g)];
    u0 = cfg.equilibrium.u0(g,:)';

    fX = @(x) room_rhs(x,u0,g,out,cfg);
    fU = @(u) room_rhs(x0,u,g,out,cfg);
    hX = @(x) sensor_output(x,cfg);

    numericA(:,:,g) = central_jacobian(fX,x0,stepX);
    numericB(:,:,g) = central_jacobian(fU,u0,stepU);
    numericCout(:,:,g) = central_jacobian(hX,x0,stepX);

    relErrorA(g) = relative_matrix_error( ...
        cont.A(:,:,g),numericA(:,:,g));
    relErrorB(g) = relative_matrix_error( ...
        cont.B(:,:,g),numericB(:,:,g));
    relErrorC(g) = relative_matrix_error( ...
        cont.Cout(:,:,g),numericCout(:,:,g));

    eigAg = eig(cont.A(:,:,g));
    maxRealEigA(g) = max(real(eigAg));
    spectralAbscissa(g) = max(abs(eigAg));
    conditionM(g) = cond(cont.M(:,:,g));

    % 令归一化状态扰动方向相同，再转换回各阶段的物理状态扰动。
    dxDirection = cont.M(:,:,g)\zDirection;
    f0 = room_rhs(x0,u0,g,out,cfg);

    for j = 1:nAlpha
        dx = alpha(j)*dxDirection;
        du = alpha(j)*uDirection;

        nonlinearIncrement = ...
            room_rhs(x0+dx,u0+du,g,out,cfg)-f0;
        linearIncrement = cont.A(:,:,g)*dx+cont.B(:,:,g)*du;

        % 在归一化输出变化率空间中度量线性化余项。
        normalizedRemainder = cont.M(:,:,g) ...
            *(nonlinearIncrement-linearIncrement);
        localRemainder(g,j) = norm(normalizedRemainder,2);
    end

    for j = 1:(nAlpha-1)
        observedOrder(g,j) = log( ...
            localRemainder(g,j)/localRemainder(g,j+1)) ...
            /log(alpha(j)/alpha(j+1));
    end
end

minimumObservedOrder = min(observedOrder,[],2);

finiteFlag = false(nStage,1);
matrixFlag = false(nStage,1);
stableFlag = false(nStage,1);
transformFlag = false(nStage,1);
quadraticFlag = false(nStage,1);

for g = 1:nStage
    finiteFlag(g) = all(isfinite(cont.A(:,:,g)),'all') ...
        && all(isfinite(cont.B(:,:,g)),'all') ...
        && all(isfinite(cont.Cout(:,:,g)),'all') ...
        && all(isfinite(cont.M(:,:,g)),'all');

    matrixFlag(g) = relErrorA(g) <= tolRelA ...
        && relErrorB(g) <= tolRelB ...
        && relErrorC(g) <= tolRelC;

    stableFlag(g) = maxRealEigA(g) < 0;
    transformFlag(g) = rank(cont.M(:,:,g)) == 3 ...
        && conditionM(g) <= tolConditionM;
    quadraticFlag(g) = minimumObservedOrder(g) >= tolOrder;
end

passFlag = finiteFlag & matrixFlag & stableFlag ...
    & transformFlag & quadraticFlag;

%% 4. 结果汇总表
summaryTable = table( ...
    cfg.stage.names,relErrorA,relErrorB,relErrorC, ...
    maxRealEigA,spectralAbscissa,conditionM,minimumObservedOrder, ...
    finiteFlag,matrixFlag,stableFlag,transformFlag,quadraticFlag,passFlag, ...
    'VariableNames',{ ...
    'Stage','RelativeErrorA','RelativeErrorB','RelativeErrorCout', ...
    'MaximumRealPartEigA','SpectralRadiusA','ConditionNumberM', ...
    'MinimumObservedOrder','FiniteMatrices','JacobianMatch', ...
    'ContinuousAStable','TransformValid','SecondOrderRemainder','Pass'});

remainderTable = table(alpha,'VariableNames',{'PerturbationScale'});
for g = 1:nStage
    remainderTable.(cfg.stage.names{g}) = localRemainder(g,:)';
end

%% 5. 保存结果
resultDir = fullfile(codeDir,'results','paper_rerun', ...
    'E2_continuous_linearization');

if exist(resultDir,'dir') ~= 7
    mkdir(resultDir);
end

writetable(summaryTable, ...
    fullfile(resultDir,'E2_linearization_validation_summary.csv'));
writetable(cont.summaryTable, ...
    fullfile(resultDir,'E2_continuous_matrix_coefficients.csv'));
writetable(remainderTable, ...
    fullfile(resultDir,'E2_local_remainder_convergence.csv'));

save(fullfile(resultDir,'E2_continuous_linearization.mat'), ...
    'cfg','cont','numericA','numericB','numericCout', ...
    'alpha','localRemainder','observedOrder','summaryTable', ...
    'remainderTable');

%% 6. 绘制局部线性化余项收敛图
fig = figure('Color','w','Position',[100 100 980 620]);
hold on;
grid on;

lineColor = lines(nStage);
for g = 1:nStage
    loglog(alpha,localRemainder(g,:),'-o', ...
        'LineWidth',1.6,'MarkerSize',6,'Color',lineColor(g,:));
end

reference = localRemainder(1,1)*(alpha/alpha(1)).^2;
loglog(alpha,reference,'k--','LineWidth',1.4);

set(gca,'XDir','reverse');
xlabel('扰动尺度 α');
ylabel('归一化局部线性化余项');
title('E2 四阶段连续局部线性化二阶余项核验');
legend([cfg.stage.names;{'二阶参考线'}],'Location','best');

saveas(fig,fullfile(resultDir,'E2_local_linearity_convergence.png'));

%% 7. 命令行输出和自动结论
fprintf('\n============================================================\n');
fprintf('E2 G1-G4四阶段连续局部线性化核验完成\n');
fprintf('参数来源：mushroom_system_parameters.m\n');
fprintf('矩阵来源：build_continuous_matrices.m\n');
fprintf('结果目录：%s\n',resultDir);
fprintf('============================================================\n');

disp(summaryTable);

fprintf('\n检查说明：\n');
fprintf('1. RelativeErrorA/B/Cout：解析矩阵与数值Jacobian误差；\n');
fprintf('2. MaximumRealPartEigA：连续A矩阵特征值最大实部，应小于0；\n');
fprintf('3. MinimumObservedOrder：局部余项收敛阶，应不小于%.1f；\n', ...
    tolOrder);
fprintf('4. Pass：以上检查同时满足时为true。\n\n');

if all(passFlag)
    fprintf('结论：PASS，四阶段连续局部线性化核验通过。\n');
    fprintf('下一步可进行E3：300 s零阶保持离散化核验。\n\n');
else
    failedStages = strjoin(cfg.stage.names(~passFlag),', ');
    error('结论：FAIL，未通过阶段：%s。',failedStages);
end
end

function J = central_jacobian(fun,x0,step)
%CENTRAL_JACOBIAN 使用中心差分计算Jacobian。
f0 = fun(x0);
J = zeros(numel(f0),numel(x0));

for j = 1:numel(x0)
    delta = zeros(size(x0));
    delta(j) = step(j);
    J(:,j) = (fun(x0+delta)-fun(x0-delta))/(2*step(j));
end
end

function value = relative_matrix_error(reference,numerical)
%RELATIVE_MATRIX_ERROR Frobenius范数相对误差。
value = norm(reference-numerical,'fro') ...
    /max(norm(reference,'fro'),eps);
end

function dx = room_rhs(x,u,g,out,cfg)
%ROOM_RHS 单仓完整非线性热湿与CO2守恒模型。
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
    cfg.physical.rhoAir_kg_m3*F_m3_s ...
      *(out.w_kg_kgda-w_kg_kgda) ...
    + mWater_kg_s + cfg.bio.mWater_kg_s(g)) ...
    /cfg.stage.Mda_kg(g);

dcCO2_kg_m3_s = ( ...
    F_m3_s*(out.cCO2_kg_m3-cCO2_kg_m3) ...
    + cfg.bio.mCO2_kg_s(g)) / cfg.physical.V_m3;

dx = [dT_C_s;dw_kg_kgda_s;dcCO2_kg_m3_s];
end

function y = sensor_output(x,cfg)
%SENSOR_OUTPUT 物理状态[T,w,c]到[T,w,CO2 ppm]的非线性输出。
T_C = x(1);
w_kg_kgda = x(2);
cCO2_kg_m3 = x(3);

CO2_ppm = cCO2_kg_m3*cfg.physical.R_J_molK*(T_C+273.15) ...
    /(cfg.physical.PIndoorRef_Pa*cfg.physical.MCO2_kg_mol)*1e6;

y = [T_C;w_kg_kgda;CO2_ppm];
end
