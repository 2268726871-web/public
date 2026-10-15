function periodic = build_periodic_cooperative_closed_loop()
%BUILD_PERIODIC_COOPERATIVE_CLOSED_LOOP
% 建立六仓周期协同闭环模型并验证局部离散稳定性
%
% 六仓阶段分配：
%   [G1,G1,G2,G3,G3,G4]
%
% 通信拓扑：
%   1-2-3-4-5-6-1无向环形拓扑
%
% 协同增益：
%   gamma = 0.35
%
% 归一化协同误差：
%   eta(k) = Sgamma*Z(k)
%
% 周期协同控制律：
%   U(k) = mathcalH*Sgamma*Z(k)
%
% 归一化闭环模型：
%   Z(k+1) = PsiGamma*Z(k)
%
% 物理坐标闭环模型：
%   X(k+1) = PhiGamma*X(k)

clc;
format long g;

%% 1. 读取六仓18维组合模型

currentFile = mfilename('fullpath');
codeDirectory = fileparts(currentFile);

resultDirectory = fullfile( ...
    codeDirectory, ...
    'results', ...
    'full_precision_matrices');

combinedResultFile = fullfile( ...
    resultDirectory, ...
    'six_room_combined_model.mat');

assert(isfile(combinedResultFile), ...
    ['没有找到六仓组合模型文件：\n', ...
     combinedResultFile, '\n', ...
     '请先运行build_six_room_combined_model。']);

data = load(combinedResultFile);

assert(isfield(data, 'combined'), ...
    '结果文件中缺少combined结构体。');

assert(isfield(data, 'cfg'), ...
    '结果文件中缺少cfg结构体。');

combined = data.combined;
cfg = data.cfg;

%% 2. 读取六仓组合矩阵

mathcalAd = combined.mathcalAd;
mathcalBd = combined.mathcalBd;
mathcalM = combined.mathcalM;
mathcalH = combined.mathcalH;

mathcalMInverse = ...
    combined.mathcalMInverse;

Abar_d = combined.Abar_d;
Bbar_d = combined.Bbar_d;

numberOfRooms = 6;
statesPerRoom = 3;
totalDimension = 18;

%% 3. 读取并核对协同增益

assert(isfield(cfg, 'control'), ...
    '唯一参数入口中缺少control字段。');

assert(isfield(cfg.control, 'gamma'), ...
    '唯一参数入口中缺少cfg.control.gamma。');

gamma = cfg.control.gamma;

assert(abs(gamma - 0.35) < 1e-12, ...
    '当前协同增益不是最终确定的0.30。');

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓周期协同闭环计算\n');
fprintf('协同增益 gamma = %.15g\n', gamma);
fprintf('=============================================\n');

%% 4. 建立六仓无向环形邻接矩阵

adjacencyMatrix = [ ...
    0 1 0 0 0 1;
    1 0 1 0 0 0;
    0 1 0 1 0 0;
    0 0 1 0 1 0;
    0 0 0 1 0 1;
    1 0 0 0 1 0];

% 每个仓的度
nodeDegree = sum(adjacencyMatrix, 2);

% 度矩阵
degreeMatrix = diag(nodeDegree);

% 拉普拉斯矩阵
laplacianMatrix = ...
    degreeMatrix - adjacencyMatrix;

fprintf('\n无向环形邻接矩阵A_ring：\n');
disp(adjacencyMatrix);

fprintf('度矩阵D：\n');
disp(degreeMatrix);

fprintf('拉普拉斯矩阵L：\n');
disp(laplacianMatrix);

%% 5. 检查无向环形拓扑

topologySymmetryError = ...
    norm(adjacencyMatrix - ...
    adjacencyMatrix.', 'fro');

numberOfEdges = ...
    sum(adjacencyMatrix(:)) / 2;

laplacianEigenvalues = ...
    sort(real(eig(laplacianMatrix)));

numberOfZeroEigenvalues = ...
    sum(abs(laplacianEigenvalues) < 1e-10);

fprintf('拉普拉斯特征值：\n');
disp(laplacianEigenvalues);

assert(topologySymmetryError < 1e-12, ...
    '邻接矩阵不是对称矩阵。');

assert(all(diag(adjacencyMatrix) == 0), ...
    '邻接矩阵对角线上出现自连接。');

assert(all(nodeDegree == 2), ...
    '无向环形拓扑中每个仓应有两个邻居。');

assert(numberOfEdges == 6, ...
    '六仓环形拓扑的边数应为6。');

assert(numberOfZeroEigenvalues == 1, ...
    '环形通信图不是连通图。');

%% 6. 建立18维协同缩放矩阵

identityRoomState = eye(statesPerRoom);

laplacianExpansion = ...
    kron(laplacianMatrix, ...
    identityRoomState);

Sgamma = ...
    eye(totalDimension) + ...
    gamma * laplacianExpansion;

SgammaEigenvalues = ...
    sort(real(eig(Sgamma)));

fprintf('\n协同缩放矩阵Sgamma的特征值：\n');
disp(SgammaEigenvalues);

assert(isequal(size(Sgamma), [18, 18]), ...
    'Sgamma维数不是18 x 18。');

assert(norm(Sgamma - Sgamma.', 'fro') < 1e-12, ...
    'Sgamma不是对称矩阵。');

assert(min(SgammaEigenvalues) > 0, ...
    'Sgamma不是正定矩阵。');

%% 7. 建立无协同独立控制闭环

% 当gamma=0时：
%
%   U(k) = mathcalH*Z(k)
%
% 物理坐标闭环：
%
%   X(k+1)
%   = [mathcalAd + mathcalBd*mathcalH*mathcalM]X(k)

PhiIndependent = ...
    mathcalAd + ...
    mathcalBd * mathcalH * mathcalM;

% 归一化坐标闭环
PsiIndependent = ...
    Abar_d + ...
    Bbar_d * mathcalH;

independentPoles = ...
    eig(PsiIndependent);

independentSpectralRadius = ...
    max(abs(independentPoles));

%% 8. 建立周期协同闭环

% 周期通信条件下，每300 s所有仓发送最新状态。
%
% 协同控制律：
%
%   U(k)
%   = mathcalH*Sgamma*Z(k)
%
%   Z(k)
%   = mathcalM*X(k)
%
% 物理坐标闭环：
%
%   X(k+1)
%   = PhiGamma*X(k)

PhiGamma = ...
    mathcalAd + ...
    mathcalBd * mathcalH * ...
    Sgamma * mathcalM;

% 归一化坐标闭环：
%
%   Z(k+1)
%   = PsiGamma*Z(k)

PsiGamma = ...
    Abar_d + ...
    Bbar_d * mathcalH * Sgamma;

%% 9. 验证物理坐标和归一化坐标一致

% 理论上应满足：
%
%   PsiGamma
%   = mathcalM*PhiGamma*mathcalM^(-1)

PsiGammaFromPhysical = ...
    mathcalM * PhiGamma * ...
    mathcalMInverse;

similarityTransformationError = ...
    norm( ...
    PsiGamma - ...
    PsiGammaFromPhysical, ...
    'fro') / ...
    (1 + norm(PsiGamma, 'fro'));

%% 10. 计算周期协同闭环特征值

periodicPoles = eig(PsiGamma);

periodicSpectralRadius = ...
    max(abs(periodicPoles));

maximumPoleMagnitude = ...
    max(abs(periodicPoles));

minimumPoleMagnitude = ...
    min(abs(periodicPoles));

numberOfUnstablePoles = ...
    sum(abs(periodicPoles) >= 1);

fprintf('\n');
fprintf('无协同独立控制谱半径 = %.15g\n', ...
    independentSpectralRadius);

fprintf('周期协同闭环谱半径   = %.15g\n', ...
    periodicSpectralRadius);

fprintf('物理与归一化闭环一致性误差 = %.15g\n', ...
    similarityTransformationError);

fprintf('模长不小于1的闭环极点个数 = %d\n', ...
    numberOfUnstablePoles);

%% 11. 输出周期协同闭环特征值

[~, poleOrder] = sort( ...
    abs(periodicPoles), ...
    'descend');

sortedPeriodicPoles = ...
    periodicPoles(poleOrder);

fprintf('\n按模长从大到小排列的周期协同闭环极点：\n');
disp(sortedPeriodicPoles);

%% 12. 自动检查

assert(all(isfinite(Sgamma(:))), ...
    'Sgamma中出现NaN或Inf。');

assert(all(isfinite(PhiGamma(:))), ...
    'PhiGamma中出现NaN或Inf。');

assert(all(isfinite(PsiGamma(:))), ...
    'PsiGamma中出现NaN或Inf。');

assert(independentSpectralRadius < 1, ...
    '无协同独立控制闭环谱半径不小于1。');

assert(similarityTransformationError < 1e-10, ...
    '物理坐标和归一化坐标闭环不一致。');

assert(periodicSpectralRadius < 1, ...
    ['gamma=0.30时周期协同闭环谱半径不小于1，', ...
     '不能判定闭环局部稳定。']);

assert(numberOfUnstablePoles == 0, ...
    '周期协同闭环中存在模长不小于1的极点。');

%% 13. 建立结果汇总表

Gamma = gamma;

IndependentSpectralRadius = ...
    independentSpectralRadius;

PeriodicCooperativeSpectralRadius = ...
    periodicSpectralRadius;

MaximumPoleMagnitude = ...
    maximumPoleMagnitude;

MinimumPoleMagnitude = ...
    minimumPoleMagnitude;

UnstablePoleCount = ...
    numberOfUnstablePoles;

SimilarityTransformationError = ...
    similarityTransformationError;

AlgebraicConnectivity = ...
    laplacianEigenvalues(2);

summaryTable = table( ...
    Gamma, ...
    IndependentSpectralRadius, ...
    PeriodicCooperativeSpectralRadius, ...
    MaximumPoleMagnitude, ...
    MinimumPoleMagnitude, ...
    UnstablePoleCount, ...
    SimilarityTransformationError, ...
    AlgebraicConnectivity);

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓周期协同闭环汇总\n');
fprintf('=============================================\n');
disp(summaryTable);

%% 14. 输出检查结果

fprintf('\n检查结果：\n');
fprintf('1. 六仓通信拓扑为连通无向环形图；\n');
fprintf('2. 每个仓均有两个通信邻居；\n');
fprintf('3. 协同增益gamma为0.30；\n');
fprintf('4. 协同缩放矩阵Sgamma为18 x 18正定矩阵；\n');
fprintf('5. 物理坐标与归一化坐标闭环一致；\n');
fprintf('6. 无协同独立控制闭环谱半径小于1；\n');
fprintf('7. 周期协同闭环谱半径小于1；\n');
fprintf('8. 六仓周期协同闭环局部离散稳定性验证通过。\n');

%% 15. 保存输出结构体

periodic.gamma = gamma;

periodic.adjacencyMatrix = ...
    adjacencyMatrix;

periodic.degreeMatrix = ...
    degreeMatrix;

periodic.laplacianMatrix = ...
    laplacianMatrix;

periodic.laplacianEigenvalues = ...
    laplacianEigenvalues;

periodic.laplacianExpansion = ...
    laplacianExpansion;

periodic.Sgamma = Sgamma;
periodic.SgammaEigenvalues = ...
    SgammaEigenvalues;

periodic.PhiIndependent = ...
    PhiIndependent;

periodic.PsiIndependent = ...
    PsiIndependent;

periodic.independentPoles = ...
    independentPoles;

periodic.independentSpectralRadius = ...
    independentSpectralRadius;

periodic.PhiGamma = PhiGamma;
periodic.PsiGamma = PsiGamma;

periodic.periodicPoles = ...
    periodicPoles;

periodic.sortedPeriodicPoles = ...
    sortedPeriodicPoles;

periodic.periodicSpectralRadius = ...
    periodicSpectralRadius;

periodic.numberOfUnstablePoles = ...
    numberOfUnstablePoles;

periodic.similarityTransformationError = ...
    similarityTransformationError;

periodic.summaryTable = summaryTable;

periodic.sampleTime_s = ...
    combined.sampleTime_s;

%% 16. 保存计算结果

save(fullfile( ...
    resultDirectory, ...
    'periodic_cooperative_closed_loop.mat'), ...
    'periodic', ...
    'combined', ...
    'cfg');

writetable(summaryTable, ...
    fullfile( ...
    resultDirectory, ...
    'periodic_cooperative_closed_loop_summary.csv'));

writematrix(adjacencyMatrix, ...
    fullfile( ...
    resultDirectory, ...
    'ring_adjacency_matrix.csv'));

writematrix(laplacianMatrix, ...
    fullfile( ...
    resultDirectory, ...
    'ring_laplacian_matrix.csv'));

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓周期协同闭环计算完成\n');
fprintf('结果目录：%s\n', resultDirectory);
fprintf('=============================================\n');

end