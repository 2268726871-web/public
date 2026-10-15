function gain = build_feedback_gains(cfg,cont)
%BUILD_FEEDBACK_GAINS 计算G1-G4全精度反馈增益矩阵。
%
% 输入：
%   cfg   唯一参数入口
%   cont  build_continuous_matrices输出
%
% 输出：
%   gain.H(:,:,g)             四阶段反馈增益矩阵
%   gain.AclLocal(:,:,g)      不考虑协同时的局部闭环矩阵
%   gain.closedLoopPoles(:,g) 局部闭环特征值
%   gain.summaryTable         增益与极点检查表
%
% 使用方法：
%   cfg = mushroom_system_parameters();
%   cont = build_continuous_matrices(cfg);
%   gain = build_feedback_gains(cfg,cont);

if nargin < 1
    cfg = mushroom_system_parameters();
end

if nargin < 2
    cont = build_continuous_matrices(cfg);
end

format long g;

nStage = cfg.stage.nStage;

H = zeros(3,3,nStage);
AclLocal = zeros(3,3,nStage);
closedLoopPoles = zeros(3,nStage);

hT = zeros(nStage,1);
hW = zeros(nStage,1);
hC = zeros(nStage,1);
hTC = zeros(nStage,1);
hWC = zeros(nStage,1);
poleError = zeros(nStage,1);

%% 1. 读取反馈设计时间常数
tauT = cfg.control.tauT_s;
tauW = cfg.control.tauW_s;
tauC = cfg.control.tauCO2_s;

desiredPoles = [ ...
    -1/tauT;
    -1/tauW;
    -1/tauC];

%% 2. 公共热容量
CT = cfg.physical.CT_J_K;

%% 3. 逐阶段计算反馈增益
for g = 1:nStage

    A_g = cont.A(:,:,g);
    B_g = cont.B(:,:,g);
    M_g = cont.M(:,:,g);

    %% 3.1 提取连续模型系数
    aT = A_g(1,1);
    aW = A_g(2,2);
    aC = A_g(3,3);

    bTF = B_g(1,3);
    bWF = B_g(2,3);
    bCF = B_g(3,3);

    %% 3.2 提取阶段归一化尺度
    sT = cfg.normalization.scaleT_C(g);
    sW = cfg.normalization.scaleW_kg_kgda(g);
    sC = cfg.normalization.scaleCO2_ppm(g);

    Mda = cfg.stage.Mda_kg(g);

    % CO2质量浓度到ppm的局部换算系数
    kPpmC = cont.Cout(3,3,g);

    %% 3.3 温度反馈增益
    hT(g) = ...
        CT*sT*(-1/tauT-aT);

    %% 3.4 含湿量反馈增益
    hW(g) = ...
        Mda*sW*(-1/tauW-aW);

    %% 3.5 CO2通风反馈增益
    assert(abs(bCF)>0, ...
        '阶段%s的通风-CO2系数为0，无法设计hC。', ...
        cfg.stage.names{g});

    hC(g) = ...
        (sC/kPpmC)*(-1/tauC-aC)/bCF;

    %% 3.6 通风对温度和含湿量的交叉补偿
    hTC(g) = -CT*bTF*hC(g);
    hWC(g) = -Mda*bWF*hC(g);

    %% 3.7 组装阶段反馈增益矩阵
    H(:,:,g) = [ ...
        hT(g), 0,     hTC(g);
        0,     hW(g), hWC(g);
        0,     0,     hC(g)];

    %% 3.8 局部闭环矩阵
    % Δxdot = AΔx + BΔu
    % Δu = H z = H M Δx
    % 所以Acl = A + BHM

    AclLocal(:,:,g) = ...
        A_g + B_g*H(:,:,g)*M_g;

    closedLoopPoles(:,g) = ...
        eig(AclLocal(:,:,g));

    %% 3.9 检查闭环极点是否等于期望极点
    actualSorted = sort(real(closedLoopPoles(:,g)));
    desiredSorted = sort(real(desiredPoles));

    poleError(g) = ...
        max(abs(actualSorted-desiredSorted));

    assert(all(isfinite(H(:,:,g)),'all'), ...
        '阶段%s的H矩阵包含无效数值。', ...
        cfg.stage.names{g});

    assert(all(real(closedLoopPoles(:,g))<0), ...
        '阶段%s的局部闭环矩阵不稳定。', ...
        cfg.stage.names{g});
end

%% 4. 增益符号检查
if any(hT>=0)
    warning('至少一个阶段的温度反馈增益hT不是负值。');
end

if any(hW>=0)
    warning('至少一个阶段的含湿量反馈增益hW不是负值。');
end

if any(hC<=0)
    warning('至少一个阶段的CO2通风增益hC不是正值。');
end

%% 5. 建立结果汇总表
stage = cfg.stage.names;

pole1 = closedLoopPoles(1,:)';
pole2 = closedLoopPoles(2,:)';
pole3 = closedLoopPoles(3,:)';

summaryTable = table( ...
    stage,hT,hTC,hW,hWC,hC, ...
    pole1,pole2,pole3,poleError, ...
    'VariableNames',{ ...
    'Stage', ...
    'hT', ...
    'hTC', ...
    'hW', ...
    'hWC', ...
    'hC', ...
    'ClosedLoopPole1', ...
    'ClosedLoopPole2', ...
    'ClosedLoopPole3', ...
    'MaxPoleError'});

%% 6. 保存输出结构
gain.H = H;
gain.AclLocal = AclLocal;
gain.closedLoopPoles = closedLoopPoles;
gain.desiredPoles = desiredPoles;
gain.poleError = poleError;
gain.summaryTable = summaryTable;

%% 7. 保存结果文件
thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile( ...
    thisDir,'results','full_precision_matrices');

if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

save(fullfile(resultDir,'feedback_gains.mat'), ...
    'gain','cfg');

writetable(summaryTable,fullfile( ...
    resultDir,'feedback_gain_summary.csv'));

%% 8. 命令行输出
fprintf('\n====================================================\n');
fprintf('G1-G4全精度反馈增益计算完成\n');
fprintf('结果目录：%s\n',resultDir);
fprintf('====================================================\n');

fprintf('\n反馈设计时间常数：\n');
fprintf('tauT   = %.15g s\n',tauT);
fprintf('tauW   = %.15g s\n',tauW);
fprintf('tauCO2 = %.15g s\n',tauC);

fprintf('\n期望连续闭环极点：\n');
disp(desiredPoles);

disp(summaryTable);

for g = 1:nStage
    fprintf('\n---------------- %s ----------------\n', ...
        cfg.stage.names{g});

    fprintf('H_%s =\n',cfg.stage.names{g});
    disp(H(:,:,g));

    fprintf('AclLocal_%s =\n',cfg.stage.names{g});
    disp(AclLocal(:,:,g));

    fprintf('eig(AclLocal_%s) =\n', ...
        cfg.stage.names{g});
    disp(closedLoopPoles(:,g));

    fprintf('最大极点误差 = %.15g\n',poleError(g));
end

fprintf('\n检查结果：\n');

if all(hT<0)
    fprintf('1. 四阶段温度反馈增益hT均小于0；\n');
else
    fprintf('1. 警告：存在非负温度反馈增益；\n');
end

if all(hW<0)
    fprintf('2. 四阶段含湿量反馈增益hW均小于0；\n');
else
    fprintf('2. 警告：存在非负含湿量反馈增益；\n');
end

if all(hC>0)
    fprintf('3. 四阶段CO2通风增益hC均大于0；\n');
else
    fprintf('3. 警告：存在非正CO2通风增益；\n');
end

fprintf('4. 四阶段局部闭环特征值实部均小于0；\n');
fprintf('5. 最大极点匹配误差 = %.15g\n',max(poleError));
fprintf('\n');

end