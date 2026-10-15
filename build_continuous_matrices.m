function cont = build_continuous_matrices(cfg)
%BUILD_CONTINUOUS_MATRICES 计算G1-G4全精度连续矩阵。
%
% 输出：
%   cont.A(:,:,g)       连续状态矩阵
%   cont.B(:,:,g)       连续输入矩阵
%   cont.Cout(:,:,g)    物理状态到传感器输出的局部转换矩阵
%   cont.Dnorm(:,:,g)   归一化尺度矩阵
%   cont.M(:,:,g)       状态归一化矩阵
%   cont.eigA(:,g)      A_g特征值
%   cont.summaryTable   四阶段关键系数汇总
%
% 使用方法：
%   cfg = mushroom_system_parameters();
%   cont = build_continuous_matrices(cfg);

if nargin < 1
    cfg = mushroom_system_parameters();
end

clc;
format long g;

nStage = cfg.stage.nStage;
nState = 3;
nInput = 3;

A = zeros(nState,nState,nStage);
B = zeros(nState,nInput,nStage);
Cout = zeros(nState,nState,nStage);
Dnorm = zeros(nState,nState,nStage);
M = zeros(nState,nState,nStage);

eigA = zeros(nState,nStage);
conditionM = zeros(nStage,1);

kPpmT = zeros(nStage,1);
kPpmC = zeros(nStage,1);

%% 1. 读取公共物理参数
rhoAir = cfg.physical.rhoAir_kg_m3;
cpAir = cfg.physical.cpAir_J_kgK;
UA = cfg.physical.UA_W_K;
CT = cfg.physical.CT_J_K;
V = cfg.physical.V_m3;

R = cfg.physical.R_J_molK;
MCO2 = cfg.physical.MCO2_kg_mol;
PIndoorRef = cfg.physical.PIndoorRef_Pa;

Tout = cfg.outdoor.nominal.T_C;
wout = cfg.outdoor.nominal.w_kg_kgda;
cout = cfg.outdoor.nominal.cCO2_kg_m3;

%% 2. 逐阶段计算全精度矩阵
for g = 1:nStage

    % 阶段目标
    Tg = cfg.stage.TTarget_C(g);
    wg = cfg.stage.wTarget_kg_kgda(g);
    cg = cfg.stage.cCO2Target_kg_m3(g);
    CppmG = cfg.stage.CO2Target_ppm(g);
    MdaG = cfg.stage.Mda_kg(g);

    % 全精度平衡通风量
    Fg = cfg.equilibrium.F0_m3_s(g);

    %% 2.1 连续状态矩阵A_g
    aT = -(UA + rhoAir*cpAir*Fg)/CT;
    aw = -(rhoAir*Fg)/MdaG;
    ac = -Fg/V;

    A(:,:,g) = [ ...
        aT, 0,  0;
        0,  aw, 0;
        0,  0,  ac];

    %% 2.2 连续输入矩阵B_g
    bTQ = 1/CT;
    bTF = rhoAir*cpAir*(Tout-Tg)/CT;

    bWm = 1/MdaG;
    bWF = rhoAir*(wout-wg)/MdaG;

    bCF = (cout-cg)/V;

    B(:,:,g) = [ ...
        bTQ, 0,   bTF;
        0,    bWm,bWF;
        0,    0,   bCF];

    %% 2.3 CO2质量浓度到ppm的局部转换
    kPpmT(g) = CppmG/(Tg+273.15);

    kPpmC(g) = ...
        R*(Tg+273.15)/(PIndoorRef*MCO2)*1e6;

    Cout(:,:,g) = [ ...
        1,        0, 0;
        0,        1, 0;
        kPpmT(g), 0, kPpmC(g)];

    %% 2.4 归一化尺度矩阵
    sT = cfg.normalization.scaleT_C(g);
    sW = cfg.normalization.scaleW_kg_kgda(g);
    sC = cfg.normalization.scaleCO2_ppm(g);

    Dnorm(:,:,g) = diag([sT,sW,sC]);

    %% 2.5 归一化矩阵M_g
    % 等价于inv(Dnorm)*Cout，但左除更加稳定。
    M(:,:,g) = Dnorm(:,:,g)\Cout(:,:,g);

    %% 2.6 数值检查
    eigA(:,g) = eig(A(:,:,g));
    conditionM(g) = cond(M(:,:,g));

    assert(all(isfinite(A(:,:,g)),'all'), ...
        '阶段%s的A矩阵包含无效数值。',cfg.stage.names{g});

    assert(all(isfinite(B(:,:,g)),'all'), ...
        '阶段%s的B矩阵包含无效数值。',cfg.stage.names{g});

    assert(all(isfinite(M(:,:,g)),'all'), ...
        '阶段%s的M矩阵包含无效数值。',cfg.stage.names{g});

    assert(rank(M(:,:,g))==3, ...
        '阶段%s的M矩阵不可逆。',cfg.stage.names{g});
end

%% 3. 检查连续A矩阵特征值
if any(real(eigA(:)) >= 0)
    warning(['至少一个连续A矩阵存在非负实部特征值，' ...
             '请检查物理参数和线性化公式。']);
end

%% 4. 建立结果汇总表
stage = cfg.stage.names;

aT = reshape(A(1,1,:),[],1);
aW = reshape(A(2,2,:),[],1);
aC = reshape(A(3,3,:),[],1);

bTQ = reshape(B(1,1,:),[],1);
bTF = reshape(B(1,3,:),[],1);
bWm = reshape(B(2,2,:),[],1);
bWF = reshape(B(2,3,:),[],1);
bCF = reshape(B(3,3,:),[],1);

summaryTable = table( ...
    stage, ...
    cfg.equilibrium.F0_m3_s, ...
    aT,aW,aC, ...
    bTQ,bTF,bWm,bWF,bCF, ...
    kPpmT,kPpmC, ...
    cfg.normalization.scaleT_C, ...
    cfg.normalization.scaleW_kg_kgda, ...
    cfg.normalization.scaleCO2_ppm, ...
    conditionM, ...
    'VariableNames',{ ...
    'Stage', ...
    'F0_m3_s', ...
    'aT_1_s','aW_1_s','aC_1_s', ...
    'bTQ','bTF','bWm','bWF','bCF', ...
    'kPpmT_ppm_C','kPpmC_ppm_per_kg_m3', ...
    'ScaleT_C','ScaleW_kg_kgda','ScaleCO2_ppm', ...
    'ConditionNumber_M'});

%% 5. 保存输出结构
cont.A = A;
cont.B = B;
cont.Cout = Cout;
cont.Dnorm = Dnorm;
cont.M = M;
cont.eigA = eigA;
cont.conditionM = conditionM;
cont.summaryTable = summaryTable;

%% 6. 输出并保存结果
thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile( ...
    thisDir,'results','full_precision_matrices');

if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

save(fullfile(resultDir, ...
    'continuous_matrices.mat'), ...
    'cont','cfg');

writetable(summaryTable,fullfile(resultDir, ...
    'continuous_matrix_summary.csv'));

fprintf('\n====================================================\n');
fprintf('G1-G4全精度连续矩阵计算完成\n');
fprintf('结果目录：%s\n',resultDir);
fprintf('====================================================\n');

disp(summaryTable);

for g = 1:nStage
    fprintf('\n---------------- %s ----------------\n', ...
        cfg.stage.names{g});

    fprintf('A_%s =\n',cfg.stage.names{g});
    disp(A(:,:,g));

    fprintf('B_%s =\n',cfg.stage.names{g});
    disp(B(:,:,g));

    fprintf('Cout_%s =\n',cfg.stage.names{g});
    disp(Cout(:,:,g));

    fprintf('M_%s =\n',cfg.stage.names{g});
    disp(M(:,:,g));

    fprintf('eig(A_%s) =\n',cfg.stage.names{g});
    disp(eigA(:,g));
end

fprintf('\n检查结果：\n');
fprintf('1. 四阶段矩阵均为有限数；\n');
fprintf('2. 四阶段M矩阵均可逆；\n');

if all(real(eigA(:))<0)
    fprintf('3. 四阶段连续A矩阵特征值实部均小于0。\n');
else
    fprintf('3. 警告：存在非负实部特征值。\n');
end

fprintf('\n');
end