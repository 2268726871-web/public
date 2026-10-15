function experiment3_zoh_discretization_validation()
%EXPERIMENT3_ZOH_DISCRETIZATION_VALIDATION
% 论文仿真实验E3：G1-G4四阶段300 s零阶保持离散化核验。
%
% 本实验验证：
%   1. Ad、Bd与增广矩阵指数给出的精确ZOH结果一致；
%   2. 离散极点与连续极点的指数映射一致；
%   3. 恒定控制输入下，离散模型的一步预测与连续模型积分一致；
%   4. 四阶段离散开环谱半径均小于1；
%   5. 四阶段局部数字反馈闭环谱半径均小于1。
%
% 说明：本实验只验证单仓局部离散模型与局部反馈稳定性，
%       不验证六仓协同控制、事件触发或设备约束。

clc;
close all;
format long g;

%% 1. 定位统一参数文件和矩阵生成函数
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

requiredFiles = { ...
    'build_continuous_matrices', ...
    'build_feedback_gains', ...
    'build_discrete_matrices'};

for iFile = 1:numel(requiredFiles)
    assert(exist(requiredFiles{iFile},'file') == 2, ...
        '未找到%s.m。',requiredFiles{iFile});
end

cfg = mushroom_system_parameters();
continuousLog = evalc('continuous = build_continuous_matrices(cfg);');
gainLog = evalc('gain = build_feedback_gains(cfg,continuous);');
discreteLog = evalc('discrete = build_discrete_matrices();');

Ts_s = cfg.time.controlPeriod_s;
nStage = cfg.stage.nStage;
nState = 3;
nInput = 3;

assert(nStage == 4,'本实验要求G1-G4四个阶段。');
assert(abs(Ts_s-300) < 1e-12, ...
    '本实验要求控制周期为300 s，当前为%.15g s。',Ts_s);
assert(abs(discrete.sampleTime_s-Ts_s) < 1e-12, ...
    '离散矩阵的采样周期与共享参数不一致。');

%% 2. 核验容差与统一测试偏差
tolAd = 1e-12;
tolBd = 1e-12;
tolPoleMapping = 1e-12;
tolOneStep = 1e-9;
tolStoredClosedLoop = 1e-12;

% 在归一化状态中设置统一初始偏差。
zInitial = [0.25;-0.20;0.15];

% 300 s内保持不变的输入偏差：[W；kg/s；m3/s]。
uHeld = [120;1.2e-5;8e-4];

%% 3. 四阶段独立核验
relativeErrorAd = zeros(nStage,1);
relativeErrorBd = zeros(nStage,1);
poleMappingError = zeros(nStage,1);
oneStepRelativeError = zeros(nStage,1);
storedClosedLoopError = zeros(nStage,1);
openLoopSpectralRadius = zeros(nStage,1);
closedLoopSpectralRadius = zeros(nStage,1);
maximumClosedPoleImaginaryPart = zeros(nStage,1);

finiteFlag = false(nStage,1);
exactZohFlag = false(nStage,1);
poleMappingFlag = false(nStage,1);
oneStepFlag = false(nStage,1);
openStableFlag = false(nStage,1);
closedStableFlag = false(nStage,1);
storedClosedLoopFlag = false(nStage,1);

odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-14);

for g = 1:nStage
    A = continuous.A(:,:,g);
    B = continuous.B(:,:,g);
    M = continuous.M(:,:,g);
    H = gain.H(:,:,g);

    Ad = discrete.Ad(:,:,g);
    Bd = discrete.Bd(:,:,g);

    % 独立参考一：直接计算exp(A*Ts)。
    AdReference = expm(A*Ts_s);

    % 独立参考二：利用A可逆时的精确ZOH等价公式。
    assert(rcond(A) > 1e-14, ...
        '阶段%s的A矩阵接近奇异，无法使用Bd等价公式。', ...
        cfg.stage.names{g});
    BdReference = A\((AdReference-eye(nState))*B);

    relativeErrorAd(g) = relative_matrix_error(Ad,AdReference);
    relativeErrorBd(g) = relative_matrix_error(Bd,BdReference);

    % 连续极点lambda应映射为离散极点exp(lambda*Ts)。
    mappedPoles = sort(real(exp(eig(A)*Ts_s)));
    actualOpenPoles = sort(real(eig(Ad)));
    poleMappingError(g) = max(abs(mappedPoles-actualOpenPoles));

    % 恒定输入下，使用ode45独立积分连续模型300 s。
    xInitial = M\zInitial;
    rhs = @(~,x) A*x+B*uHeld;
    [~,xTrajectory] = ode45(rhs,[0 Ts_s],xInitial,odeOptions);
    xContinuousEnd = xTrajectory(end,:)';
    xDiscreteEnd = Ad*xInitial+Bd*uHeld;

    normalizedDifference = M*(xContinuousEnd-xDiscreteEnd);
    normalizedReference = M*xContinuousEnd;
    oneStepRelativeError(g) = norm(normalizedDifference,2) ...
        /max(norm(normalizedReference,2),1e-12);

    % 300 s采样、周期内保持输入不变时的局部数字反馈闭环。
    AclDigital = Ad+Bd*H*M;
    storedClosedLoopError(g) = relative_matrix_error( ...
        discrete.AclDiscrete(:,:,g),AclDigital);

    openPoles = eig(Ad);
    closedPoles = eig(AclDigital);
    openLoopSpectralRadius(g) = max(abs(openPoles));
    closedLoopSpectralRadius(g) = max(abs(closedPoles));
    maximumClosedPoleImaginaryPart(g) = ...
        max(abs(imag(closedPoles)));

    finiteFlag(g) = all(isfinite(Ad),'all') ...
        && all(isfinite(Bd),'all') ...
        && all(isfinite(AclDigital),'all');

    exactZohFlag(g) = relativeErrorAd(g) <= tolAd ...
        && relativeErrorBd(g) <= tolBd;
    poleMappingFlag(g) = poleMappingError(g) <= tolPoleMapping;
    oneStepFlag(g) = oneStepRelativeError(g) <= tolOneStep;
    openStableFlag(g) = openLoopSpectralRadius(g) < 1;
    closedStableFlag(g) = closedLoopSpectralRadius(g) < 1;
    storedClosedLoopFlag(g) = ...
        storedClosedLoopError(g) <= tolStoredClosedLoop;
end

passFlag = finiteFlag & exactZohFlag & poleMappingFlag ...
    & oneStepFlag & openStableFlag & closedStableFlag ...
    & storedClosedLoopFlag;

%% 4. 结果汇总表
summaryTable = table( ...
    cfg.stage.names,relativeErrorAd,relativeErrorBd,poleMappingError, ...
    oneStepRelativeError,openLoopSpectralRadius, ...
    closedLoopSpectralRadius,maximumClosedPoleImaginaryPart, ...
    storedClosedLoopError,finiteFlag,exactZohFlag,poleMappingFlag, ...
    oneStepFlag,openStableFlag,closedStableFlag, ...
    storedClosedLoopFlag,passFlag, ...
    'VariableNames',{ ...
    'Stage','RelativeErrorAd','RelativeErrorBd','PoleMappingError', ...
    'OneStepRelativeError','OpenLoopSpectralRadius', ...
    'ClosedLoopSpectralRadius','MaximumClosedPoleImaginaryPart', ...
    'StoredClosedLoopError','FiniteMatrices','ExactZOH', ...
    'PoleMappingMatch','OneStepMatch','OpenLoopStable', ...
    'LocalDigitalClosedLoopStable','StoredClosedLoopMatch','Pass'});

%% 5. 保存正式重跑结果
resultDir = fullfile(codeDir,'results','paper_rerun', ...
    'E3_zoh_discretization');

if exist(resultDir,'dir') ~= 7
    mkdir(resultDir);
end

writetable(summaryTable, ...
    fullfile(resultDir,'E3_zoh_validation_summary.csv'));
writetable(discrete.summaryTable, ...
    fullfile(resultDir,'E3_discrete_matrix_summary.csv'));

save(fullfile(resultDir,'E3_zoh_discretization.mat'), ...
    'cfg','continuous','gain','discrete','summaryTable', ...
    'zInitial','uHeld','continuousLog','gainLog','discreteLog');

%% 6. 绘制谱半径与一步预测误差
fig = figure('Color','w','Position',[100 100 1080 720]);

subplot(2,1,1);
bar([openLoopSpectralRadius,closedLoopSpectralRadius]);
hold on;
yline(1,'r--','稳定边界','LineWidth',1.4);
grid on;
set(gca,'XTick',1:nStage,'XTickLabel',cfg.stage.names);
ylabel('谱半径');
title('E3 四阶段300 s离散模型稳定性核验');
legend({'离散开环','局部数字反馈闭环'},'Location','best');
ylim([0 1.05]);

subplot(2,1,2);
% 对数坐标不能显示0；仅绘图时将0和更小数值放在10^-16显示下限。
plotOneStepError = max(oneStepRelativeError,1e-16);
semilogy(1:nStage,plotOneStepError,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;
yline(tolOneStep,'r--','通过上限','LineWidth',1.4);
grid on;
set(gca,'XTick',1:nStage,'XTickLabel',cfg.stage.names);
xlabel('生长阶段');
ylabel('300 s一步预测相对误差');
title('恒定输入下连续积分与离散一步预测的一致性');
ylim([1e-16 1e-8]);

saveas(fig,fullfile(resultDir,'E3_zoh_validation.png'));

%% 7. 命令行结论
fprintf('\n============================================================\n');
fprintf('E3 G1-G4四阶段300 s零阶保持离散化核验完成\n');
fprintf('采样周期：%.15g s\n',Ts_s);
fprintf('参数来源：mushroom_system_parameters.m\n');
fprintf('结果目录：%s\n',resultDir);
fprintf('============================================================\n');

disp(summaryTable);

fprintf('\n检查说明：\n');
fprintf('1. ExactZOH：Ad、Bd与独立精确公式一致；\n');
fprintf('2. PoleMappingMatch：连续极点与离散极点指数映射一致；\n');
fprintf('3. OneStepMatch：恒定输入下300 s连续积分与离散预测一致；\n');
fprintf('4. OpenLoopStable：离散开环谱半径小于1；\n');
fprintf('5. LocalDigitalClosedLoopStable：局部数字反馈闭环谱半径小于1；\n');
fprintf('6. Pass：以上检查同时满足时为true。\n\n');

if all(passFlag)
    fprintf('结论：PASS，四阶段300 s零阶保持离散化核验通过。\n');
    fprintf('离散模型可用于后续六仓控制与通信实验。\n\n');
else
    failedStages = strjoin(cfg.stage.names(~passFlag),', ');
    error('结论：FAIL，未通过阶段：%s。',failedStages);
end
end

function value = relative_matrix_error(reference,numerical)
%RELATIVE_MATRIX_ERROR Frobenius范数相对误差。
value = norm(reference-numerical,'fro') ...
    /max(norm(reference,'fro'),eps);
end
