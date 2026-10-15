function combined = build_six_room_combined_model()
%BUILD_SIX_ROOM_COMBINED_MODEL
% 建立六个蘑菇仓的18维离散组合模型
%
% 六仓阶段分配：
%   [G1,G1,G2,G3,G3,G4]
%
% 单仓状态：
%   Delta_x_i = [Delta_T_i; Delta_w_i; Delta_c_i]
%
% 单仓输入：
%   Delta_u_i = [Delta_Q_i; Delta_mw_i; Delta_F_i]
%
% 六仓组合模型：
%   X(k+1) = mathcal_Ad*X(k) + mathcal_Bd*U(k)
%
% 归一化状态：
%   Z(k) = mathcal_M*X(k)
%
% 本函数只建立18维组合模型。
% 暂不加入环形拓扑和协同增益gamma。

clc;
format long g;

%% 1. 定位前面保存的全精度离散结果

currentFile = mfilename('fullpath');
codeDirectory = fileparts(currentFile);

resultDirectory = fullfile( ...
    codeDirectory, ...
    'results', ...
    'full_precision_matrices');

discreteResultFile = fullfile( ...
    resultDirectory, ...
    'full_precision_discrete_matrices.mat');

assert(isfile(discreteResultFile), ...
    ['没有找到全精度离散结果文件：\n', ...
     discreteResultFile, '\n', ...
     '请先运行 build_discrete_matrices。']);

%% 2. 读取经过验证的数据

data = load(discreteResultFile);

assert(isfield(data, 'cfg'), ...
    '离散结果文件中缺少cfg。');

assert(isfield(data, 'continuous'), ...
    '离散结果文件中缺少continuous。');

assert(isfield(data, 'gain'), ...
    '离散结果文件中缺少gain。');

assert(isfield(data, 'discrete'), ...
    '离散结果文件中缺少discrete。');

cfg = data.cfg;
continuous = data.continuous;
gain = data.gain;
discrete = data.discrete;

%% 3. 六仓阶段分配

% 阶段编号：
%   1代表G1
%   2代表G2
%   3代表G3
%   4代表G4

stageIndex = [1, 1, 2, 3, 3, 4];

expectedStageIndex = [1, 1, 2, 3, 3, 4];

assert(isequal(stageIndex, expectedStageIndex), ...
    '六仓阶段分配不正确。');

numberOfRooms = 6;
statesPerRoom = 3;
inputsPerRoom = 3;

totalNumberOfStates = ...
    numberOfRooms * statesPerRoom;

totalNumberOfInputs = ...
    numberOfRooms * inputsPerRoom;

assert(totalNumberOfStates == 18, ...
    '六仓组合状态维数不是18。');

assert(totalNumberOfInputs == 18, ...
    '六仓组合输入维数不是18。');

%% 4. 初始化六仓组合矩阵

mathcalAd = zeros( ...
    totalNumberOfStates, ...
    totalNumberOfStates);

mathcalBd = zeros( ...
    totalNumberOfStates, ...
    totalNumberOfInputs);

mathcalM = zeros( ...
    totalNumberOfStates, ...
    totalNumberOfStates);

mathcalH = zeros( ...
    totalNumberOfInputs, ...
    totalNumberOfStates);

%% 5. 根据阶段分配逐仓组装矩阵

for room = 1:numberOfRooms

    stage = stageIndex(room);

    stateIndex = ...
        (room - 1) * statesPerRoom + ...
        (1:statesPerRoom);

    inputIndex = ...
        (room - 1) * inputsPerRoom + ...
        (1:inputsPerRoom);

    % 离散状态矩阵
    mathcalAd(stateIndex, stateIndex) = ...
        discrete.Ad(:, :, stage);

    % 离散输入矩阵
    mathcalBd(stateIndex, inputIndex) = ...
        discrete.Bd(:, :, stage);

    % 归一化和传感器转换矩阵
    mathcalM(stateIndex, stateIndex) = ...
        continuous.M(:, :, stage);

    % 反馈增益矩阵
    mathcalH(inputIndex, stateIndex) = ...
        gain.H(:, :, stage);

end

%% 6. 建立归一化坐标下的六仓模型

% 定义：
%
%   Z(k) = mathcalM*X(k)
%
% 因此：
%
%   X(k) = mathcalM^(-1)*Z(k)
%
% 代入物理坐标模型：
%
%   Z(k+1)
%   = mathcalM*mathcalAd*mathcalM^(-1)*Z(k)
%   + mathcalM*mathcalBd*U(k)
%
% 定义：
%
%   Abar_d = mathcalM*mathcalAd*mathcalM^(-1)
%   Bbar_d = mathcalM*mathcalBd

assert(rcond(mathcalM) > 1e-14, ...
    '六仓组合归一化矩阵mathcalM不可逆或病态。');

mathcalMInverse = ...
    mathcalM \ eye(totalNumberOfStates);

Abar_d = ...
    mathcalM * mathcalAd * mathcalMInverse;

Bbar_d = ...
    mathcalM * mathcalBd;

%% 7. 生成六仓阶段分配表

RoomNumber = (1:numberOfRooms).';

StageIndex = stageIndex(:);

StageName = strings(numberOfRooms, 1);

for room = 1:numberOfRooms
    StageName(room) = ...
        string(cfg.stage.names{stageIndex(room)});
end

roomStageTable = table( ...
    RoomNumber, ...
    StageIndex, ...
    StageName, ...
    'VariableNames', { ...
    'RoomNumber', ...
    'StageIndex', ...
    'StageName'});

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓阶段分配\n');
fprintf('=============================================\n');
disp(roomStageTable);

%% 8. 组合矩阵维数检查

fprintf('\n');
fprintf('六仓组合矩阵维数：\n');

fprintf('size(mathcalAd) = %d x %d\n', ...
    size(mathcalAd, 1), ...
    size(mathcalAd, 2));

fprintf('size(mathcalBd) = %d x %d\n', ...
    size(mathcalBd, 1), ...
    size(mathcalBd, 2));

fprintf('size(mathcalM)  = %d x %d\n', ...
    size(mathcalM, 1), ...
    size(mathcalM, 2));

fprintf('size(mathcalH)  = %d x %d\n', ...
    size(mathcalH, 1), ...
    size(mathcalH, 2));

fprintf('size(Abar_d)    = %d x %d\n', ...
    size(Abar_d, 1), ...
    size(Abar_d, 2));

fprintf('size(Bbar_d)    = %d x %d\n', ...
    size(Bbar_d, 1), ...
    size(Bbar_d, 2));

%% 9. 建立参考块对角矩阵进行核对

AdCell = cell(numberOfRooms, 1);
BdCell = cell(numberOfRooms, 1);
MCell = cell(numberOfRooms, 1);
HCell = cell(numberOfRooms, 1);

for room = 1:numberOfRooms

    stage = stageIndex(room);

    AdCell{room} = ...
        discrete.Ad(:, :, stage);

    BdCell{room} = ...
        discrete.Bd(:, :, stage);

    MCell{room} = ...
        continuous.M(:, :, stage);

    HCell{room} = ...
        gain.H(:, :, stage);

end

mathcalAdReference = blkdiag(AdCell{:});
mathcalBdReference = blkdiag(BdCell{:});
mathcalMReference = blkdiag(MCell{:});
mathcalHReference = blkdiag(HCell{:});

AdAssemblyError = ...
    norm(mathcalAd - ...
    mathcalAdReference, 'fro');

BdAssemblyError = ...
    norm(mathcalBd - ...
    mathcalBdReference, 'fro');

MAssemblyError = ...
    norm(mathcalM - ...
    mathcalMReference, 'fro');

HAssemblyError = ...
    norm(mathcalH - ...
    mathcalHReference, 'fro');

%% 10. 归一化坐标变换一致性检查

% 构造确定性的测试状态和测试输入
testX = (1:totalNumberOfStates).' / ...
    totalNumberOfStates;

testU = (totalNumberOfInputs:-1:1).' / ...
    totalNumberOfInputs;

% 先在物理坐标中计算下一时刻，再转换为Z
nextZFromPhysicalModel = ...
    mathcalM * ...
    (mathcalAd * testX + ...
     mathcalBd * testU);

% 直接在归一化坐标中计算下一时刻
testZ = mathcalM * testX;

nextZFromNormalizedModel = ...
    Abar_d * testZ + ...
    Bbar_d * testU;

coordinateTransformationError = ...
    norm( ...
    nextZFromPhysicalModel - ...
    nextZFromNormalizedModel, 2) / ...
    (1 + norm(nextZFromPhysicalModel, 2));

%% 11. 六仓开环谱半径

combinedOpenLoopPoles = eig(mathcalAd);

combinedOpenLoopSpectralRadius = ...
    max(abs(combinedOpenLoopPoles));

%% 12. 自动检查

assert(isequal(size(mathcalAd), [18, 18]), ...
    'mathcalAd维数不是18 x 18。');

assert(isequal(size(mathcalBd), [18, 18]), ...
    'mathcalBd维数不是18 x 18。');

assert(isequal(size(mathcalM), [18, 18]), ...
    'mathcalM维数不是18 x 18。');

assert(isequal(size(mathcalH), [18, 18]), ...
    'mathcalH维数不是18 x 18。');

assert(all(isfinite(mathcalAd(:))), ...
    'mathcalAd中存在NaN或Inf。');

assert(all(isfinite(mathcalBd(:))), ...
    'mathcalBd中存在NaN或Inf。');

assert(all(isfinite(mathcalM(:))), ...
    'mathcalM中存在NaN或Inf。');

assert(all(isfinite(mathcalH(:))), ...
    'mathcalH中存在NaN或Inf。');

assert(AdAssemblyError < 1e-12, ...
    'mathcalAd块对角组装错误。');

assert(BdAssemblyError < 1e-12, ...
    'mathcalBd块对角组装错误。');

assert(MAssemblyError < 1e-12, ...
    'mathcalM块对角组装错误。');

assert(HAssemblyError < 1e-12, ...
    'mathcalH块对角组装错误。');

assert(coordinateTransformationError < 1e-10, ...
    '物理坐标与归一化坐标模型不一致。');

assert(combinedOpenLoopSpectralRadius < 1, ...
    '六仓组合开环谱半径不小于1。');

%% 13. 输出检查结果

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓18维组合模型检查结果\n');
fprintf('=============================================\n');

fprintf('mathcalAd组装误差 = %.15g\n', ...
    AdAssemblyError);

fprintf('mathcalBd组装误差 = %.15g\n', ...
    BdAssemblyError);

fprintf('mathcalM组装误差  = %.15g\n', ...
    MAssemblyError);

fprintf('mathcalH组装误差  = %.15g\n', ...
    HAssemblyError);

fprintf('坐标变换一致性误差 = %.15g\n', ...
    coordinateTransformationError);

fprintf('六仓开环谱半径 = %.15g\n', ...
    combinedOpenLoopSpectralRadius);

fprintf('\n检查结果：\n');
fprintf('1. 六仓阶段分配为[G1,G1,G2,G3,G3,G4]；\n');
fprintf('2. 四个组合矩阵均为18 x 18；\n');
fprintf('3. 六仓矩阵均按照阶段正确组装；\n');
fprintf('4. 六仓组合归一化矩阵可逆；\n');
fprintf('5. 物理坐标与归一化坐标模型一致；\n');
fprintf('6. 六仓18维组合模型建立通过。\n');

%% 14. 保存输出结构体

combined.stageIndex = stageIndex;
combined.roomStageTable = roomStageTable;

combined.mathcalAd = mathcalAd;
combined.mathcalBd = mathcalBd;
combined.mathcalM = mathcalM;
combined.mathcalH = mathcalH;

combined.mathcalMInverse = ...
    mathcalMInverse;

combined.Abar_d = Abar_d;
combined.Bbar_d = Bbar_d;

combined.openLoopPoles = ...
    combinedOpenLoopPoles;

combined.openLoopSpectralRadius = ...
    combinedOpenLoopSpectralRadius;

combined.AdAssemblyError = ...
    AdAssemblyError;

combined.BdAssemblyError = ...
    BdAssemblyError;

combined.MAssemblyError = ...
    MAssemblyError;

combined.HAssemblyError = ...
    HAssemblyError;

combined.coordinateTransformationError = ...
    coordinateTransformationError;

combined.stateDimension = ...
    totalNumberOfStates;

combined.inputDimension = ...
    totalNumberOfInputs;

combined.sampleTime_s = ...
    discrete.sampleTime_s;

%% 15. 保存结果

save(fullfile( ...
    resultDirectory, ...
    'six_room_combined_model.mat'), ...
    'combined', ...
    'cfg');

writetable(roomStageTable, ...
    fullfile( ...
    resultDirectory, ...
    'six_room_stage_assignment.csv'));

fprintf('\n');
fprintf('=============================================\n');
fprintf('六仓18维组合模型计算完成\n');
fprintf('结果目录：%s\n', resultDirectory);
fprintf('=============================================\n');

end