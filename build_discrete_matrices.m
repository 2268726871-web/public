function discrete = build_discrete_matrices()
%BUILD_DISCRETE_MATRICES
% G1-G4连续模型的300 s零阶保持精确离散化
%
% 连续模型：
%   Delta_x_dot(t) = A_g*Delta_x(t) + B_g*Delta_u(t)
%
% 离散模型：
%   Delta_x(k+1) = Ad_g*Delta_x(k) + Bd_g*Delta_u(k)
%
% 控制输入在每个300 s控制周期内保持不变。
%
% 输出结构体：
%   discrete.Ad
%   discrete.Bd
%   discrete.AclDiscrete
%   discrete.openLoopPoles
%   discrete.closedLoopPoles
%   discrete.spectralRadiusOpen
%   discrete.spectralRadiusClosed
%   discrete.summaryTable

clc;
format long g;

%% 1. 读取唯一参数、连续模型和反馈增益

cfg = mushroom_system_parameters();

continuous = build_continuous_matrices();

gain = build_feedback_gains();

Ts = cfg.time.controlPeriod_s;

assert(abs(Ts - 300) < 1e-12, ...
    '控制周期不是300 s，请检查唯一参数入口。');

fprintf('\n');
fprintf('============================================================\n');
fprintf('开始进行G1-G4零阶保持精确离散化\n');
fprintf('控制采样周期 Ts = %.15g s\n', Ts);
fprintf('============================================================\n');

%% 2. 矩阵维数

numberOfStages = 4;
numberOfStates = 3;
numberOfInputs = 3;

% 四阶段离散矩阵
Ad = zeros(numberOfStates, ...
           numberOfStates, ...
           numberOfStages);

Bd = zeros(numberOfStates, ...
           numberOfInputs, ...
           numberOfStages);

% 四阶段离散闭环矩阵
AclDiscrete = zeros(numberOfStates, ...
                    numberOfStates, ...
                    numberOfStages);

% 特征值
openLoopPoles = zeros(numberOfStates, ...
                      numberOfStages);

closedLoopPoles = zeros(numberOfStates, ...
                        numberOfStages);

% 谱半径
spectralRadiusOpen = zeros(numberOfStages, 1);
spectralRadiusClosed = zeros(numberOfStages, 1);

% 数值检查
zohCheckError = zeros(numberOfStages, 1);
maximumPoleImaginaryPart = zeros(numberOfStages, 1);

%% 3. 四阶段零阶保持精确离散化

for g = 1:numberOfStages

    % 当前阶段连续矩阵
    Ag = continuous.A(:, :, g);
    Bg = continuous.B(:, :, g);

    % 当前阶段归一化矩阵和反馈增益
    Mg = continuous.M(:, :, g);
    Hg = gain.H(:, :, g);

    %----------------------------------------------------------
    % 构造增广矩阵：
    %
    %             [Ag  Bg]
    % Augmented = [      ]
    %             [0   0 ]
    %
    % 计算：
    %
    % expm(Augmented*Ts) = [Ad  Bd]
    %                      [0    I]
    %
    % 由此一次性获得精确Ad和Bd。
    %----------------------------------------------------------

    augmentedMatrix = ...
        [Ag, Bg;
         zeros(numberOfInputs, ...
         numberOfStates + numberOfInputs)];

    augmentedExponential = ...
        expm(augmentedMatrix * Ts);

    Ad(:, :, g) = ...
        augmentedExponential( ...
        1:numberOfStates, ...
        1:numberOfStates);

    Bd(:, :, g) = ...
        augmentedExponential( ...
        1:numberOfStates, ...
        numberOfStates + 1 : ...
        numberOfStates + numberOfInputs);

    %----------------------------------------------------------
    % 离散闭环模型
    %
    % 反馈控制律：
    %
    % Delta_u(k) = Hg*Mg*Delta_x(k)
    %
    % 代入离散模型：
    %
    % Delta_x(k+1)
    % = [Ad + Bd*Hg*Mg]Delta_x(k)
    %----------------------------------------------------------

    AclDiscrete(:, :, g) = ...
        Ad(:, :, g) + ...
        Bd(:, :, g) * Hg * Mg;

    % 开环离散特征值
    openLoopPoles(:, g) = ...
        eig(Ad(:, :, g));

    % 闭环离散特征值
    closedLoopPoles(:, g) = ...
        eig(AclDiscrete(:, :, g));

    % 谱半径：所有特征值模长的最大值
    spectralRadiusOpen(g) = ...
        max(abs(openLoopPoles(:, g)));

    spectralRadiusClosed(g) = ...
        max(abs(closedLoopPoles(:, g)));

    % 检查闭环极点是否存在明显虚部
    maximumPoleImaginaryPart(g) = ...
        max(abs(imag(closedLoopPoles(:, g))));

    %----------------------------------------------------------
    % 使用等价公式核对Bd
    %
    % 如果Ag可逆：
    %
    % Bd = Ag^(-1)*(Ad-I)*Bg
    %
    % MATLAB中使用Ag\，不直接计算inv(Ag)。
    %----------------------------------------------------------

    if rcond(Ag) > 1e-14

        BdEquivalent = ...
            Ag \ ((Ad(:, :, g) - ...
            eye(numberOfStates)) * Bg);

        zohCheckError(g) = ...
            norm(Bd(:, :, g) - ...
            BdEquivalent, 'fro');

    else

        zohCheckError(g) = NaN;

    end

end

%% 4. 连续期望极点对应的参考离散极点

% 连续期望极点：
%   -1/tauT
%   -1/tauW
%   -1/tauCO2
%
% 相应的参考离散极点：
%   exp(lambdaContinuous*Ts)
%
% 注意：
% 这些只是参考值。
% 实际控制指令每300 s更新一次，并在周期内保持不变，
% 因此实际采样闭环极点不要求与参考值完全相等。

referenceDiscretePoles = ...
    exp(gain.desiredPoles * Ts);

fprintf('\n');
fprintf('连续期望极点对应的参考离散极点：\n');
disp(referenceDiscretePoles);

%% 5. 输出四阶段计算结果

for g = 1:numberOfStages

    stageName = cfg.stage.names{g};

    fprintf('\n');
    fprintf('-------------------- %s --------------------\n', ...
        stageName);

    fprintf('Ad_%s =\n', stageName);
    disp(Ad(:, :, g));

    fprintf('Bd_%s =\n', stageName);
    disp(Bd(:, :, g));

    fprintf('AclDiscrete_%s =\n', stageName);
    disp(AclDiscrete(:, :, g));

    fprintf('eig(Ad_%s) =\n', stageName);
    disp(openLoopPoles(:, g));

    fprintf('eig(AclDiscrete_%s) =\n', stageName);
    disp(closedLoopPoles(:, g));

    fprintf('开环谱半径 = %.15g\n', ...
        spectralRadiusOpen(g));

    fprintf('闭环谱半径 = %.15g\n', ...
        spectralRadiusClosed(g));

    fprintf('Bd等价公式检查误差 = %.15g\n', ...
        zohCheckError(g));

end

%% 6. 建立结果汇总表

Stage = cfg.stage.names(:);

summaryTable = table( ...
    Stage, ...
    spectralRadiusOpen, ...
    spectralRadiusClosed, ...
    zohCheckError, ...
    maximumPoleImaginaryPart, ...
    'VariableNames', { ...
    'Stage', ...
    'OpenLoopSpectralRadius', ...
    'ClosedLoopSpectralRadius', ...
    'BdEquivalentFormulaError', ...
    'MaximumPoleImaginaryPart'});

fprintf('\n');
fprintf('================ 离散化结果汇总 ================\n');
disp(summaryTable);

%% 7. 自动检查

assert(all(isfinite(Ad(:))), ...
    'Ad中出现NaN或Inf。');

assert(all(isfinite(Bd(:))), ...
    'Bd中出现NaN或Inf。');

assert(all(isfinite(AclDiscrete(:))), ...
    '离散闭环矩阵中出现NaN或Inf。');

assert(all(spectralRadiusOpen < 1), ...
    '至少一个阶段的离散开环谱半径不小于1。');

assert(all(spectralRadiusClosed < 1), ...
    '至少一个阶段的离散闭环谱半径不小于1。');

assert(all(maximumPoleImaginaryPart < 1e-10), ...
    '离散闭环极点出现了不可忽略的虚部。');

validZohError = ...
    zohCheckError(~isnan(zohCheckError));

assert(all(validZohError < 1e-10), ...
    'Bd的两种计算方法不一致，请检查离散化过程。');

fprintf('\n');
fprintf('检查结果：\n');
fprintf('1. 四阶段Ad和Bd均为有限数；\n');
fprintf('2. 四阶段离散开环谱半径均小于1；\n');
fprintf('3. 四阶段离散闭环谱半径均小于1；\n');
fprintf('4. Bd两种计算方式一致；\n');
fprintf('5. 300 s零阶保持精确离散化通过。\n');

%% 8. 保存输出结构体

discrete.Ad = Ad;
discrete.Bd = Bd;

discrete.AclDiscrete = AclDiscrete;

discrete.openLoopPoles = openLoopPoles;
discrete.closedLoopPoles = closedLoopPoles;

discrete.spectralRadiusOpen = ...
    spectralRadiusOpen;

discrete.spectralRadiusClosed = ...
    spectralRadiusClosed;

discrete.zohCheckError = ...
    zohCheckError;

discrete.maximumPoleImaginaryPart = ...
    maximumPoleImaginaryPart;

discrete.referenceDiscretePoles = ...
    referenceDiscretePoles;

discrete.sampleTime_s = Ts;
discrete.summaryTable = summaryTable;

%% 9. 保存MAT和CSV文件

currentFile = mfilename('fullpath');
codeDirectory = fileparts(currentFile);

resultDirectory = fullfile( ...
    codeDirectory, ...
    'results', ...
    'full_precision_matrices');

if ~exist(resultDirectory, 'dir')
    mkdir(resultDirectory);
end

save(fullfile( ...
    resultDirectory, ...
    'full_precision_discrete_matrices.mat'), ...
    'discrete', ...
    'cfg', ...
    'continuous', ...
    'gain');

writetable(summaryTable, ...
    fullfile( ...
    resultDirectory, ...
    'discrete_matrix_summary.csv'));

fprintf('\n');
fprintf('============================================================\n');
fprintf('G1-G4全精度离散矩阵计算完成\n');
fprintf('结果目录：%s\n', resultDirectory);
fprintf('============================================================\n');

end