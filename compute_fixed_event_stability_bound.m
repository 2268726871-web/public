function fixedStability = compute_fixed_event_stability_bound()
%COMPUTE_FIXED_EVENT_STABILITY_BOUND
% 基于全精度六仓周期协同闭环矩阵，计算固定事件触发
% Lyapunov充分稳定上界thetaStab
%
% 事件触发闭环：
%   Z(k+1) = PsiGamma*Z(k) + Ge*E(k)
%
% 离散Lyapunov方程：
%   PsiGamma'*P*PsiGamma - P = -Q
%
% 本研究选择：
%   Q = I_18
%
% 固定触发误差约束：
%   ||E(k)|| <= theta*||Z(k)||
%
% 稳定充分条件：
%   -lambda_min(Q)
%   + 2*alphaFixed*theta
%   + betaFixed*theta^2 < 0

clc;
format long g;

%% 1. 读取周期协同闭环全精度结果

currentFile = mfilename('fullpath');
codeDirectory = fileparts(currentFile);

resultDirectory = fullfile( ...
    codeDirectory, ...
    'results', ...
    'full_precision_matrices');

periodicResultFile = fullfile( ...
    resultDirectory, ...
    'periodic_cooperative_closed_loop.mat');

assert(isfile(periodicResultFile), ...
    ['没有找到周期协同闭环结果文件：\n', ...
     periodicResultFile, '\n', ...
     '请先运行build_periodic_cooperative_closed_loop。']);

data = load(periodicResultFile);

assert(isfield(data, 'periodic'), ...
    '结果文件中缺少periodic结构体。');

assert(isfield(data, 'combined'), ...
    '结果文件中缺少combined结构体。');

assert(isfield(data, 'cfg'), ...
    '结果文件中缺少cfg结构体。');

periodic = data.periodic;
combined = data.combined;
cfg = data.cfg;

%% 2. 读取全精度闭环矩阵

PsiGamma = periodic.PsiGamma;

Bbar_d = combined.Bbar_d;
mathcalH = combined.mathcalH;
Sgamma = periodic.Sgamma;

gamma = periodic.gamma;

numberOfStates = size(PsiGamma, 1);

assert(numberOfStates == 18, ...
    '周期协同闭环矩阵不是18维。');

assert(abs(gamma - 0.35) < 1e-12, ...
    '协同增益不是最终确定的0.35。');

%% 3. 再次检查周期协同闭环稳定性

periodicPoles = eig(PsiGamma);

periodicSpectralRadius = ...
    max(abs(periodicPoles));

assert(periodicSpectralRadius < 1, ...
    ['周期协同闭环谱半径不小于1，', ...
     '不能求取正定离散Lyapunov矩阵。']);

fprintf('\n');
fprintf('====================================================\n');
fprintf('固定事件触发Lyapunov稳定上界复算\n');
fprintf('协同增益gamma = %.15g\n', gamma);
fprintf('周期协同闭环谱半径 = %.15g\n', ...
    periodicSpectralRadius);
fprintf('====================================================\n');

%% 4. 建立事件误差输入矩阵Ge

% 事件触发控制使用：
%
%   hatZ(k) = Z(k) + E(k)
%
%   U(k) = mathcalH*Sgamma*hatZ(k)
%
% 所以：
%
%   Z(k+1)
%   = PsiGamma*Z(k) + Ge*E(k)
%
% 其中：
%
%   Ge = Bbar_d*mathcalH*Sgamma

Ge = ...
    Bbar_d * mathcalH * Sgamma;

assert(isequal(size(Ge), [18, 18]), ...
    '事件误差输入矩阵Ge不是18 x 18。');

assert(all(isfinite(Ge(:))), ...
    '事件误差输入矩阵Ge中出现NaN或Inf。');

%% 5. 选择离散Lyapunov方程中的Q矩阵

% 这里选择Q = I_18。
% 这是理论分析中的标准选择，不是设备参数。

Qf = eye(numberOfStates);

qMinimum = ...
    min(real(eig(Qf)));

assert(qMinimum > 0, ...
    'Lyapunov方程中的Qf不是正定矩阵。');

%% 6. 求解离散Lyapunov方程

% 需要求解：
%
%   PsiGamma'*Pf*PsiGamma - Pf = -Qf
%
% 等价于：
%
%   Pf - PsiGamma'*Pf*PsiGamma = Qf
%
% 使用向量化公式：
%
% vec(A*P*B)
% = (B'*kron*A)*vec(P)
%
% 得到：
%
% [I - kron(PsiGamma',PsiGamma')]*vec(Pf)
% = vec(Qf)
%
% 这种方法不依赖Control System Toolbox。

lyapunovOperator = ...
    eye(numberOfStates^2) - ...
    kron(PsiGamma.', PsiGamma.');

assert(rcond(lyapunovOperator) > 1e-14, ...
    '离散Lyapunov线性方程病态或不可逆。');

Pvector = ...
    lyapunovOperator \ Qf(:);

Pf = ...
    reshape(Pvector, ...
    numberOfStates, ...
    numberOfStates);

% 消除浮点计算造成的极小非对称误差
Pf = (Pf + Pf.') / 2;

%% 7. 检查Lyapunov方程残差和Pf正定性

lyapunovResidualMatrix = ...
    PsiGamma.' * Pf * PsiGamma ...
    - Pf + Qf;

lyapunovResidual = ...
    norm(lyapunovResidualMatrix, 'fro') / ...
    (1 + norm(Qf, 'fro'));

PfSymmetryError = ...
    norm(Pf - Pf.', 'fro');

PfEigenvalues = ...
    sort(real(eig(Pf)));

minimumEigenvaluePf = ...
    min(PfEigenvalues);

maximumEigenvaluePf = ...
    max(PfEigenvalues);

conditionNumberPf = ...
    cond(Pf);

assert(lyapunovResidual < 1e-9, ...
    '离散Lyapunov方程残差过大。');

assert(PfSymmetryError < 1e-10, ...
    '计算得到的Pf不是对称矩阵。');

assert(minimumEigenvaluePf > 0, ...
    '计算得到的Pf不是正定矩阵。');

%% 8. 计算alphaFixed和betaFixed

crossMatrix = ...
    PsiGamma.' * Pf * Ge;

eventQuadraticMatrix = ...
    Ge.' * Pf * Ge;

% 对称化，用于消除浮点误差
eventQuadraticMatrix = ...
    (eventQuadraticMatrix + ...
     eventQuadraticMatrix.') / 2;

alphaFixed = ...
    norm(crossMatrix, 2);

betaFixed = ...
    norm(eventQuadraticMatrix, 2);

assert(isfinite(alphaFixed) && ...
       alphaFixed >= 0, ...
    'alphaFixed计算错误。');

assert(isfinite(betaFixed) && ...
       betaFixed >= 0, ...
    'betaFixed计算错误。');

%% 9. 计算固定事件触发稳定上界

% 一般形式：
%
% betaFixed*theta^2
% + 2*alphaFixed*theta
% - qMinimum < 0
%
% 正根即充分稳定上界。

if betaFixed > 1e-14

    % 直接公式
    thetaStabDirect = ...
        (sqrt(alphaFixed^2 + ...
        betaFixed * qMinimum) ...
        - alphaFixed) / betaFixed;

    % 数值稳定的等价公式
    thetaStab = ...
        qMinimum / ...
        (sqrt(alphaFixed^2 + ...
        betaFixed * qMinimum) ...
        + alphaFixed);

else

    assert(alphaFixed > 0, ...
        'alphaFixed和betaFixed同时接近零。');

    thetaStabDirect = ...
        qMinimum / (2 * alphaFixed);

    thetaStab = thetaStabDirect;

end

thetaFormulaConsistencyError = ...
    abs(thetaStab - thetaStabDirect);

assert(isfinite(thetaStab), ...
    'thetaStab不是有限数。');

assert(thetaStab > 0, ...
    'thetaStab不大于零。');

%% 10. 检查稳定不等式

% 在理论上界处，二次函数应接近零。

coefficientAtBoundary = ...
    -qMinimum ...
    + 2 * alphaFixed * thetaStab ...
    + betaFixed * thetaStab^2;

% 使用上界的99%进行内部检查。
% 0.99只用于验证不等式，不是最终实验阈值。

thetaInsideCheck = ...
    0.99 * thetaStab;

coefficientInsideBound = ...
    -qMinimum ...
    + 2 * alphaFixed * thetaInsideCheck ...
    + betaFixed * thetaInsideCheck^2;

assert(abs(coefficientAtBoundary) < 1e-8, ...
    'thetaStab代入二次边界后不接近零。');

assert(coefficientInsideBound < 0, ...
    '小于thetaStab的测试阈值未满足充分稳定条件。');

%% 11. 输出计算结果

fprintf('\n');
fprintf('离散Lyapunov方程结果：\n');
fprintf('Qf最小特征值 = %.15g\n', ...
    qMinimum);

fprintf('Pf最小特征值 = %.15g\n', ...
    minimumEigenvaluePf);

fprintf('Pf最大特征值 = %.15g\n', ...
    maximumEigenvaluePf);

fprintf('Pf条件数 = %.15g\n', ...
    conditionNumberPf);

fprintf('Lyapunov方程相对残差 = %.15g\n', ...
    lyapunovResidual);

fprintf('\n');
fprintf('事件误差稳定系数：\n');
fprintf('alphaFixed = %.15g\n', ...
    alphaFixed);

fprintf('betaFixed  = %.15g\n', ...
    betaFixed);

fprintf('\n');
fprintf('固定事件触发充分稳定上界：\n');
fprintf('thetaStab = %.15g\n', ...
    thetaStab);

fprintf('边界处二次函数值 = %.15g\n', ...
    coefficientAtBoundary);

fprintf('0.99*thetaStab处二次函数值 = %.15g\n', ...
    coefficientInsideBound);

%% 12. 建立结果汇总表

Gamma = gamma;

PeriodicSpectralRadius = ...
    periodicSpectralRadius;

QMinimumEigenvalue = ...
    qMinimum;

PMinimumEigenvalue = ...
    minimumEigenvaluePf;

PMaximumEigenvalue = ...
    maximumEigenvaluePf;

PConditionNumber = ...
    conditionNumberPf;

LyapunovResidual = ...
    lyapunovResidual;

AlphaFixed = ...
    alphaFixed;

BetaFixed = ...
    betaFixed;

ThetaStab = ...
    thetaStab;

CoefficientAtBoundary = ...
    coefficientAtBoundary;

summaryTable = table( ...
    Gamma, ...
    PeriodicSpectralRadius, ...
    QMinimumEigenvalue, ...
    PMinimumEigenvalue, ...
    PMaximumEigenvalue, ...
    PConditionNumber, ...
    LyapunovResidual, ...
    AlphaFixed, ...
    BetaFixed, ...
    ThetaStab, ...
    CoefficientAtBoundary);

fprintf('\n');
fprintf('====================================================\n');
fprintf('固定事件触发稳定上界汇总\n');
fprintf('====================================================\n');
disp(summaryTable);

%% 13. 输出检查结果

fprintf('\n检查结果：\n');
fprintf('1. 周期协同闭环谱半径小于1；\n');
fprintf('2. 事件误差输入矩阵Ge为18 x 18有限矩阵；\n');
fprintf('3. 离散Lyapunov方程求解残差满足要求；\n');
fprintf('4. Lyapunov矩阵Pf为对称正定矩阵；\n');
fprintf('5. alphaFixed和betaFixed均为有限非负数；\n');
fprintf('6. thetaStab为正的有限数；\n');
fprintf('7. 小于thetaStab的测试阈值满足充分稳定条件；\n');
fprintf('8. 固定事件触发Lyapunov稳定上界复算通过。\n');

%% 14. 保存输出结构体

fixedStability.gamma = gamma;

fixedStability.PsiGamma = ...
    PsiGamma;

fixedStability.Ge = Ge;

fixedStability.Qf = Qf;
fixedStability.Pf = Pf;

fixedStability.PfEigenvalues = ...
    PfEigenvalues;

fixedStability.minimumEigenvaluePf = ...
    minimumEigenvaluePf;

fixedStability.maximumEigenvaluePf = ...
    maximumEigenvaluePf;

fixedStability.conditionNumberPf = ...
    conditionNumberPf;

fixedStability.lyapunovResidual = ...
    lyapunovResidual;

fixedStability.alphaFixed = ...
    alphaFixed;

fixedStability.betaFixed = ...
    betaFixed;

fixedStability.thetaStab = ...
    thetaStab;

fixedStability.thetaFormulaConsistencyError = ...
    thetaFormulaConsistencyError;

fixedStability.coefficientAtBoundary = ...
    coefficientAtBoundary;

fixedStability.thetaInsideCheck = ...
    thetaInsideCheck;

fixedStability.coefficientInsideBound = ...
    coefficientInsideBound;

fixedStability.summaryTable = ...
    summaryTable;

%% 15. 保存结果文件

save(fullfile( ...
    resultDirectory, ...
    'fixed_event_stability_bound.mat'), ...
    'fixedStability', ...
    'periodic', ...
    'combined', ...
    'cfg');

writetable(summaryTable, ...
    fullfile( ...
    resultDirectory, ...
    'fixed_event_stability_bound_summary.csv'));

fprintf('\n');
fprintf('====================================================\n');
fprintf('固定事件触发稳定上界复算完成\n');
fprintf('理论充分稳定范围：0 < theta_f < %.15g\n', ...
    thetaStab);
fprintf('结果目录：%s\n', resultDirectory);
fprintf('====================================================\n');

end