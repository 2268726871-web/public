function experiment2_collaboration_ablation_4units_gamma035()
%EXPERIMENT2_COLLABORATION_ABLATION 实验二：六仓协同作用消融实验。
%
% 对照原则：
%   A. 独立控制：gamma = 0，各仓只使用本仓状态；
%   B. 协同控制：gamma = 0.30，六仓经无向环形网络交换归一化误差。
% 除gamma和是否需要仓间通信外，两组的物理模型、初始状态、控制器增益、
% 设备约束、积分步长、控制周期、气象扰动与仿真时长完全相同。
%
% 运行方法：
%   将MATLAB当前目录切换到本文件所在目录，执行：
%       experiment2_collaboration_ablation
%
% 结果目录：
%   本文件所在目录/results/experiment2_collaboration_ablation/
%
% 数值性质：
%   [L] 文献/数据表原始值；[D] 派生计算值；[A] 算法设计值。
%   程序会把关键数值及其性质写入parameter_provenance.csv。

clc;
close all;

%% 1. 时钟和气象文件
dt = 60;                         % s，用户指定
controlPeriod = 300;             % s，用户指定
weatherPeriod = 3600;            % s，用户指定
tEnd = 72*3600;                  % s，用户指定
t = (0:dt:tEnd)';
nStep = numel(t);
nRoom = 6;
nControlStep = controlPeriod/dt;

if mod(controlPeriod,dt) ~= 0 || mod(weatherPeriod,dt) ~= 0
    error('控制周期和气象周期必须是积分步长的整数倍。');
end

weatherFile = ...
    '/Users/bystander/Desktop/蘑菇仓/南宁2023逐时室外气象与CO2背景数据_重建版.xlsx';
if ~isfile(weatherFile)
    error(['未找到气象文件：\n%s\n请修改代码第35行附近的weatherFile。'], ...
        weatherFile);
end

%% 2. 物理模型参数（继承step1与实验一）
par.rhoAir = 1.2;                % kg/m^3
par.cpAir = 1006;                % J/(kg K)
par.UA = 31.5;                   % W/K
par.CT = 6.85e6;                 % J/K
par.V = 104;                     % m^3
par.R = 8.314462618;             % J/(mol K)
par.MCO2 = 44.01e-3;             % kg/mol
par.PIndoorRef = 101325;         % Pa

%% 3. 读取全年气象，并用可复现规则选出典型72 h
% Excel“第四部分模型输入”：
% A北京时间，B室外温度，C室外RH，D表压(kPa)，E绝对压力(Pa)，
% F 2023全球CO2背景代理，G城市CO2情景。本实验采用B/C/E/F。
raw = readmatrix(weatherFile, 'Sheet', '第四部分模型输入', ...
    'Range', 'A2:G8761');
if size(raw,1) ~= 8760 || size(raw,2) < 6
    error('气象表应含8760行且至少6列，当前尺寸为%d x %d。', ...
        size(raw,1), size(raw,2));
end
if any(~isfinite(raw(:,1:6)), 'all')
    error('气象表A:F范围含缺失或非数值项，请先检查原始工作簿。');
end

allDate = datetime(raw(:,1), 'ConvertFrom', 'excel');
allTout = raw(:,2);
allRH = raw(:,3);
allP = raw(:,5);
allCO2 = raw(:,6);
allWout = humidity_ratio_from_t_rh_p(allTout, allRH, allP);

% [A] 只比较从每日00:00开始的连续72 h候选窗口。
% 特征=[mean(T),mean(w),std(T),std(w)]。均值目标取全年8760 h均值，
% 波动目标取所有候选窗口波动的中位数；四项用候选IQR稳健标准化后
% 计算等权平方距离，距离最小者定义为“典型72 h”。这既与实验一
% 使用的全年标称均值可比，又避免按控制结果挑选天气。
[weather, selection] = select_typical_72h( ...
    allDate, allTout, allRH, allP, allCO2, allWout);

weatherHourIndex = min(floor(t/weatherPeriod)+1, 72);
outT = weather.T(weatherHourIndex);
outW = weather.w(weatherHourIndex);
outP = weather.P(weatherHourIndex);
outCO2ppm = weather.CO2ppm(weatherHourIndex);
outC = ppm_to_mass(outCO2ppm, outT, outP, par);

%% 4. 阶段目标、扰动与六仓配置
stageNames = {'G1';'G2';'G3';'G4'};
roomNames = {'1号仓';'2号仓';'3号仓';'4号仓';'5号仓';'6号仓'};
roomStage = [1;1;2;3;3;4];

T0 = [23;24;16;16];
w0 = [0.0168318;0.0179081;0.0096345;0.0096345];
Cppm0 = [12500;12500;900;900];
Mda = [120.7;120.1;125.0;125.0];

Qbio = [170;246;246;189];
mbioW = [7.87e-5;8.36e-5;1.53e-4;1.53e-4];
mbioC = [4.73e-5;6.83e-5;6.83e-5;5.25e-5];

TMin = [22;23;15;15];
TMax = [24;25;17;17];
wMin = [0.0149619;0.0159232;0.0084889;0.0084889];
wMax = [0.0179081;0.0190463;0.0108933;0.0108933];
CppmMin = [10000;10000;800;800];
CppmMax = [15000;15000;1000;1000];

scaleT = [1;1;1;1];
scaleW = [0.0010763;0.0011382;0.0011456;0.0011456];
scaleCppm = [2500;2500;100;100];

initial.T = [23.30,22.70,24.20,16.40,15.80,16.30];
% 2号仓必须为0.0165089；作废值0.0171547未使用。
initial.w = [0.0166165,0.0165089,0.0175666,0.0098636,0.0094054,0.0099782];
initial.Cppm = [13500,13000,13750,960,930,940];

%% 5. 平衡前馈、设备约束与控制器
% 下列标称值只保留作step1一致性参照。动态气象实验中，两组均在每小时
% 气象更新后，用同一组平衡方程重新派生uFF(Tout,wout,Pout,CO2out)。
% 这是已知室外扰动的模型前馈，不是新增反馈增益，也不是事件触发机制。
uFFNominal = [-147.543, -7.29243e-5, 0.002160593;
              -187.493, -7.11874e-5, 0.003130736;
              -1011.754,-5.93899e-4, 0.073932066;
              -823.821, -4.91905e-4, 0.056829187];

device.coolRated_W = 7310;
device.heatRated_W = 9810;
device.humMax_kg_h = 2.40;
device.dehPerUnit_kg_h = 0.8375;
device.nDehumidifier = 4;
device.fanLow_m3_h = 600;
device.fanHigh_m3_h = 1000;
device.coolPowerRated_kW = 2.00;
device.heatPowerRated_kW = 3.00;
device.humPowerRated_kW = 0.210;
device.dehPowerPerUnit_kW = 0.50;
device.fanLowPower_kW = 0.203;
device.fanHighPower_kW = 0.550;

limits.QMin = -device.coolRated_W;
limits.QMax = device.heatRated_W;
limits.mwMin = -(device.nDehumidifier*device.dehPerUnit_kg_h)/3600;
limits.mwMax = device.humMax_kg_h/3600;
limits.FMin = 0;
limits.FMax = device.fanHigh_m3_h/3600;

% 六仓无向环：1-2-3-4-5-6-1。
A = [0 1 0 0 0 1;
     1 0 1 0 0 0;
     0 1 0 1 0 0;
     0 0 1 0 1 0;
     0 0 0 1 0 1;
     1 0 0 0 1 0];
L = diag(sum(A,2))-A;

H = zeros(3,3,4);
H(:,:,1) = [-917.281,0,18.640;
             0,-1.0547e-4,6.2690e-5;
             0,0,0.023451];
H(:,:,2) = [-916.109,0,46.554;
             0,-1.0964e-4,9.2194e-5;
             0,0,0.023253];
H(:,:,3) = [-830.638,0,-63.979;
             0,-1.7697e-5,-4.9839e-5;
             0,0,0.008357];
H(:,:,4) = [-851.285,0,-90.268;
             0,-4.1209e-5,-7.0317e-5;
             0,0,0.011791];

gammaIndependent = 0.0;          % [A] 消融：删除协同项
gammaCollaborative = 0.35;       % [A] 与实验一一致

model = struct('par',par,'T0',T0,'w0',w0,'Cppm0',Cppm0,'Mda',Mda, ...
    'Qbio',Qbio,'mbioW',mbioW,'mbioC',mbioC,'TMin',TMin,'TMax',TMax, ...
    'wMin',wMin,'wMax',wMax,'CppmMin',CppmMin,'CppmMax',CppmMax, ...
    'scaleT',scaleT,'scaleW',scaleW,'scaleCppm',scaleCppm, ...
    'roomStage',roomStage,'uFFNominal',uFFNominal,'H',H,'L',L,'A',A, ...
    'limits',limits,'device',device);
clock = struct('dt',dt,'controlPeriod',controlPeriod, ...
    'nControlStep',nControlStep,'t',t,'nStep',nStep,'nRoom',nRoom);
outdoor = struct('T',outT,'w',outW,'P',outP,'Cppm',outCO2ppm,'c',outC);

%% 6. 两组闭环仿真
fprintf('\n正在运行A组：独立控制（gamma=0）...\n');
independent = simulate_case('Independent', gammaIndependent, false, ...
    initial, model, clock, outdoor);
fprintf('正在运行B组：协同控制（gamma=0.35）...\n');
collaborative = simulate_case('Collaborative', gammaCollaborative, true, ...
    initial, model, clock, outdoor);

%% 7. 统一指标计算
independent = analyze_case(independent, model, clock);
collaborative = analyze_case(collaborative, model, clock);

comparison = build_comparison_table(independent, collaborative);
effect = build_effect_table(independent, collaborative);
roomComparison = build_room_comparison( ...
    independent, collaborative, roomNames, stageNames, roomStage);

% 核心判据只用于自动概括，不把它冒充文献标准：
% 协同组图一致性RMSE和同阶段对RMSE均下降，且安全越界不增加。
ablationSupportsCollaboration = ...
    collaborative.metric.GraphDisagreementRMSE < ...
        independent.metric.GraphDisagreementRMSE && ...
    collaborative.metric.SameStagePairRMSE < ...
        independent.metric.SameStagePairRMSE && ...
    collaborative.metric.SafetyViolationSamples <= ...
        independent.metric.SafetyViolationSamples;

%% 8. 保存结果
thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(thisDir,'results', ...
    'experiment2_collaboration_ablation_4units_gamma035');
if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

writetable(comparison,fullfile(resultDir,'comparison_summary.csv'));
writetable(effect,fullfile(resultDir,'collaboration_effect.csv'));
writetable(roomComparison,fullfile(resultDir,'room_comparison.csv'));

selectedWeatherTable = table((0:71)',weather.Date,weather.T,weather.RH, ...
    weather.P,weather.w,weather.CO2ppm, ...
    'VariableNames',{'HourIndex','DateTime_Beijing','T_out_C','RH_out_pct', ...
    'P_out_Pa','w_out_kg_kgda','CO2_background_ppm'});
writetable(selectedWeatherTable,fullfile(resultDir,'selected_weather_72h.csv'));

selectionTable = table(selection.SelectedStart,selection.SelectedEnd, ...
    selection.CandidateCount,selection.Score,selection.MeanT, ...
    selection.MeanW,selection.StdT,selection.StdW,selection.MinT, ...
    selection.MaxT,selection.MinRH,selection.MaxRH, ...
    'VariableNames',{'SelectedStart_Beijing','SelectedEnd_Beijing', ...
    'CandidateCount','RobustDistanceScore','MeanT_C','MeanW_kg_kgda', ...
    'StdT_C','StdW_kg_kgda','MinT_C','MaxT_C','MinRH_pct','MaxRH_pct'});
writetable(selectionTable,fullfile(resultDir,'weather_selection_summary.csv'));

write_case_outputs(resultDir,'independent',independent,clock,outdoor, ...
    roomNames,stageNames,roomStage,model);
write_case_outputs(resultDir,'collaborative',collaborative,clock,outdoor, ...
    roomNames,stageNames,roomStage,model);
write_provenance(fullfile(resultDir,'parameter_provenance.csv'), ...
    weatherFile,selection,dt,controlPeriod,weatherPeriod,tEnd, ...
    gammaIndependent,gammaCollaborative,device,limits);

save(fullfile(resultDir,'experiment2_collaboration_ablation.mat'), ...
    'independent','collaborative','comparison','effect','roomComparison', ...
    'weather','selection','model','clock','outdoor', ...
    'ablationSupportsCollaboration');

try
    make_figures(resultDir,clock,weather,independent,collaborative,roomNames);
catch ME
    warning('数值结果已保存，但当前MATLAB图形环境导出PNG失败：%s',ME.message);
end

%% 9. 命令行摘要
fprintf('\n============================================================\n');
fprintf('实验二：六仓协同作用消融实验完成\n');
fprintf('唯一控制差异：gamma %.2f（独立） vs %.2f（协同）\n', ...
    gammaIndependent,gammaCollaborative);
fprintf('动态气象：%s 至 %s（北京时间，共72 h）\n', ...
    char(selection.SelectedStart),char(selection.SelectedEnd));
fprintf('候选窗口：%d个；稳健距离得分：%.8f\n', ...
    selection.CandidateCount,selection.Score);
fprintf('结果目录：%s\n',resultDir);
fprintf('============================================================\n');
disp(comparison);
disp(effect);
if ablationSupportsCollaboration
    fprintf(['结论：协同组同时降低图一致性误差和同阶段仓对误差，', ...
        '且未增加安全越界；本次消融支持保留协同项。\n\n']);
else
    fprintf(['结论：本次结果未同时满足全部协同证据判据。', ...
        '请依据collaboration_effect.csv逐项讨论，不应选择性汇报。\n\n']);
end
end

function result = simulate_case(caseName,gamma,usesCommunication, ...
    initial,model,clock,outdoor)
% 两组调用同一个仿真函数，避免对照组出现隐藏实现差异。
t = clock.t;
dt = clock.dt;
nStep = clock.nStep;
nRoom = clock.nRoom;

T = zeros(nStep,nRoom);
w = zeros(nStep,nRoom);
c = zeros(nStep,nRoom);
T(1,:) = initial.T;
w(1,:) = initial.w;
for i = 1:nRoom
    c(1,i) = ppm_to_mass(initial.Cppm(i),T(1,i), ...
        model.par.PIndoorRef,model.par);
end

Qraw = zeros(nStep,nRoom);
mwRaw = zeros(nStep,nRoom);
Fraw = zeros(nStep,nRoom);
Qcmd = zeros(nStep,nRoom);
mwCmd = zeros(nStep,nRoom);
Fcmd = zeros(nStep,nRoom);
satQ = false(nStep,nRoom);
satW = false(nStep,nRoom);
satF = false(nStep,nRoom);
satDetail = false(nStep,nRoom,6);
isControlUpdate = false(nStep,1);
isBroadcast = false(nStep,nRoom);

uRawHold = zeros(nRoom,3);
uHold = zeros(nRoom,3);
satHold = false(nRoom,6);

for k = 1:(nStep-1)
    out = struct('T',outdoor.T(k),'w',outdoor.w(k), ...
        'P',outdoor.P(k),'Cppm',outdoor.Cppm(k),'c',outdoor.c(k));
    if mod(k-1,clock.nControlStep) == 0
        % ================================================================
        % 控制律第一次实际使用：状态误差如何变成设备指令
        % ================================================================
        % 1) 本仓状态先变为无量纲误差：
        %    z_i=[(T_i-T*_g)/s_Tg,(w_i-w*_g)/s_wg,
        %         (C_i-C*_g)/s_Cg]。
        %    s取目标至最近安全边界的距离，所以温度、含湿量和CO2可以
        %    在同一协同公式中比较，而不会因单位不同失衡。
        %
        % 2) 消融变量只出现在下式：
        %    eta_i=z_i+gamma*sum_{j in N_i}(z_i-z_j)。
        %    A组gamma=0，故eta_i=z_i，各仓只纠正自身偏差，不广播；
        %    B组gamma=0.30，六仓每300 s广播最新z，第二项抑制环形邻仓
        %    归一化偏差不一致。跨阶段比较的是“相对各自目标的偏差”，
        %    不是把不同生长阶段强行控制到相同温度或CO2。
        %
        % 3) 误差转为净物理输入：
        %    u_raw,i=u_ff,g(Tout,wout,Pout,Cout)+H_g*eta_i，
        %    u=[Q_T,m_w,F_v]。u_ff由step1的三条平衡方程随逐时气象
        %    重新派生；两组使用完全相同的u_ff，因此不会污染协同消融。
        %    eta_T>0表示偏热；H的温度对角元为负，使Q_T减小，即增加
        %    制冷。eta_w>0表示偏湿，湿量对角元为负，使m_w减小，即
        %    增加除湿。eta_C>0表示CO2偏高，通风增益为正，使F_v增大。
        %    H中的交叉项补偿通风同时造成的显热和湿量耦合。
        %
        % 4) u_raw经过设备上下限后得到u_cmd。Q_T用正负号唯一选择
        %    加热/制冷，m_w用正负号唯一选择加湿/除湿，因此两组都满足
        %    制冷-加热和加湿-除湿互锁；F_v限制为0~1000 m3/h。
        [uRawHold,uHold,satHold] = controller_step( ...
            T(k,:),w(k,:),c(k,:),model,gamma,out);
        isControlUpdate(k) = true;
        if usesCommunication
            isBroadcast(k,:) = true;
        end
    end

    Qraw(k,:) = uRawHold(:,1)';
    mwRaw(k,:) = uRawHold(:,2)';
    Fraw(k,:) = uRawHold(:,3)';
    Qcmd(k,:) = uHold(:,1)';
    mwCmd(k,:) = uHold(:,2)';
    Fcmd(k,:) = uHold(:,3)';
    satQ(k,:) = (satHold(:,1)|satHold(:,2))';
    satW(k,:) = (satHold(:,3)|satHold(:,4))';
    satF(k,:) = (satHold(:,5)|satHold(:,6))';
    satDetail(k,:,:) = reshape(satHold,[1,nRoom,6]);

    for i = 1:nRoom
        g = model.roomStage(i);
        xk = [T(k,i);w(k,i);c(k,i)];
        uRoom = uHold(i,:)';
        bio = [model.Qbio(g);model.mbioW(g);model.mbioC(g)];
        k1 = room_rhs(xk,uRoom,bio,model.Mda(g),out,model.par);
        k2 = room_rhs(xk+0.5*dt*k1,uRoom,bio,model.Mda(g),out,model.par);
        k3 = room_rhs(xk+0.5*dt*k2,uRoom,bio,model.Mda(g),out,model.par);
        k4 = room_rhs(xk+dt*k3,uRoom,bio,model.Mda(g),out,model.par);
        xNext = xk+dt/6*(k1+2*k2+2*k3+k4);
        T(k+1,i) = xNext(1);
        w(k+1,i) = xNext(2);
        c(k+1,i) = max(xNext(3),0);
    end
end

% 记录72 h端点的控制计算，但不再积分。
out = struct('T',outdoor.T(end),'w',outdoor.w(end), ...
    'P',outdoor.P(end),'Cppm',outdoor.Cppm(end),'c',outdoor.c(end));
[uRawHold,uHold,satHold] = controller_step( ...
    T(end,:),w(end,:),c(end,:),model,gamma,out);
isControlUpdate(end) = true;
if usesCommunication
    isBroadcast(end,:) = true;
end
Qraw(end,:) = uRawHold(:,1)';
mwRaw(end,:) = uRawHold(:,2)';
Fraw(end,:) = uRawHold(:,3)';
Qcmd(end,:) = uHold(:,1)';
mwCmd(end,:) = uHold(:,2)';
Fcmd(end,:) = uHold(:,3)';
satQ(end,:) = (satHold(:,1)|satHold(:,2))';
satW(end,:) = (satHold(:,3)|satHold(:,4))';
satF(end,:) = (satHold(:,5)|satHold(:,6))';
satDetail(end,:,:) = reshape(satHold,[1,nRoom,6]);

result = struct('Name',caseName,'gamma',gamma, ...
    'usesCommunication',usesCommunication,'T',T,'w',w,'c',c, ...
    'Qraw',Qraw,'mwRaw',mwRaw,'Fraw',Fraw, ...
    'Qcmd',Qcmd,'mwCmd',mwCmd,'Fcmd',Fcmd, ...
    'satQ',satQ,'satW',satW,'satF',satF,'satDetail',satDetail, ...
    'isControlUpdate',isControlUpdate,'isBroadcast',isBroadcast);
end

function [uRaw,uApplied,sat] = controller_step(Trow,wrow,crow,model,gamma,out)
nRoom = numel(model.roomStage);
z = zeros(nRoom,3);
for i = 1:nRoom
    g = model.roomStage(i);
    CppmNow = mass_to_ppm(crow(i),Trow(i), ...
        model.par.PIndoorRef,model.par);
    z(i,:) = [(Trow(i)-model.T0(g))/model.scaleT(g), ...
              (wrow(i)-model.w0(g))/model.scaleW(g), ...
              (CppmNow-model.Cppm0(g))/model.scaleCppm(g)];
end

eta = z+gamma*(model.L*z);
uFFNow = equilibrium_feedforward(model,out);
uRaw = zeros(nRoom,3);
for i = 1:nRoom
    g = model.roomStage(i);
    uRaw(i,:) = uFFNow(g,:)+(model.H(:,:,g)*eta(i,:)')';
end

uApplied = uRaw;
uApplied(:,1) = min(max(uRaw(:,1),model.limits.QMin),model.limits.QMax);
uApplied(:,2) = min(max(uRaw(:,2),model.limits.mwMin),model.limits.mwMax);
uApplied(:,3) = min(max(uRaw(:,3),model.limits.FMin),model.limits.FMax);
sat = [uRaw(:,1)<model.limits.QMin,uRaw(:,1)>model.limits.QMax, ...
       uRaw(:,2)<model.limits.mwMin,uRaw(:,2)>model.limits.mwMax, ...
       uRaw(:,3)<model.limits.FMin,uRaw(:,3)>model.limits.FMax];
end

function uFF = equilibrium_feedforward(model,out)
% 直接由step1使用的平衡方程派生当前逐时室外状态下的前馈输入。
nStage = numel(model.T0);
cTarget = zeros(nStage,1);
for g = 1:nStage
    cTarget(g) = ppm_to_mass(model.Cppm0(g),model.T0(g), ...
        model.par.PIndoorRef,model.par);
end
F0 = model.mbioC./(cTarget-out.c);
mw0 = model.par.rhoAir.*F0.*(model.w0-out.w)-model.mbioW;
Q0 = (model.par.UA+model.par.rhoAir*model.par.cpAir.*F0) ...
    .*(model.T0-out.T)-model.Qbio;
uFF = [Q0,mw0,F0];
end

function dx = room_rhs(x,u,bio,Mda,out,par)
T = x(1);
w = x(2);
c = x(3);
QT = u(1);
mw = u(2);
Fv = u(3);
dT = (par.UA*(out.T-T)+par.rhoAir*par.cpAir*Fv*(out.T-T) ...
    +QT+bio(1))/par.CT;
dw = (par.rhoAir*Fv*(out.w-w)+mw+bio(2))/Mda;
dc = (Fv*(out.c-c)+bio(3))/par.V;
dx = [dT;dw;dc];
end

function result = analyze_case(result,model,clock)
nStep = clock.nStep;
nRoom = clock.nRoom;
Cppm = mass_to_ppm(result.c,result.T,model.par.PIndoorRef,model.par);

targetT = zeros(1,nRoom);
targetW = zeros(1,nRoom);
targetC = zeros(1,nRoom);
safeTMin = zeros(1,nRoom);
safeTMax = zeros(1,nRoom);
safeWMin = zeros(1,nRoom);
safeWMax = zeros(1,nRoom);
safeCMin = zeros(1,nRoom);
safeCMax = zeros(1,nRoom);
normT = zeros(1,nRoom);
normW = zeros(1,nRoom);
normC = zeros(1,nRoom);
for i = 1:nRoom
    g = model.roomStage(i);
    targetT(i)=model.T0(g); targetW(i)=model.w0(g); targetC(i)=model.Cppm0(g);
    safeTMin(i)=model.TMin(g); safeTMax(i)=model.TMax(g);
    safeWMin(i)=model.wMin(g); safeWMax(i)=model.wMax(g);
    safeCMin(i)=model.CppmMin(g); safeCMax(i)=model.CppmMax(g);
    normT(i)=model.scaleT(g); normW(i)=model.scaleW(g);
    normC(i)=model.scaleCppm(g);
end

dT = result.T-targetT;
dW = result.w-targetW;
dC = Cppm-targetC;
zT = dT./normT;
zW = dW./normW;
zC = dC./normC;
z = cat(3,zT,zW,zC);

violTLow = result.T<safeTMin;
violTHigh = result.T>safeTMax;
violWLow = result.w<safeWMin;
violWHigh = result.w>safeWMax;
violCLow = Cppm<safeCMin;
violCHigh = Cppm>safeCMax;
anyViolation = violTLow|violTHigh|violWLow|violWHigh|violCLow|violCHigh;

[~,~,~,~,~,~,~,~,~,~,~,devicePower] = decode_devices( ...
    result.Qcmd,result.mwCmd,result.Fcmd,model.device);
energyPerRoom = sum(devicePower(1:end-1,:),1)'*clock.dt/3600;

edge = [1 2;2 3;3 4;4 5;5 6;6 1];
graphDisagreement = zeros(nStep,1);
for k = 1:nStep
    accum = [];
    for e = 1:size(edge,1)
        dz = squeeze(z(k,edge(e,1),:)-z(k,edge(e,2),:));
        accum = [accum;dz(:)]; %#ok<AGROW>
    end
    graphDisagreement(k) = sqrt(mean(accum.^2));
end

sameStagePair = [1 2;4 5];
sameStageDisagreement = zeros(nStep,1);
for k = 1:nStep
    accum = [];
    for p = 1:size(sameStagePair,1)
        dz = squeeze(z(k,sameStagePair(p,1),:)- ...
            z(k,sameStagePair(p,2),:));
        accum = [accum;dz(:)]; %#ok<AGROW>
    end
    sameStageDisagreement(k) = sqrt(mean(accum.^2));
end

first6 = clock.t<=6*3600;
after6 = clock.t>=6*3600;
controlIndex = find(result.isControlUpdate);
trackingTime = sqrt(mean(reshape(z.^2,nStep,[]),2));

nBroadcast = nnz(result.isBroadcast);
nReceive = 2*nBroadcast; % 环形网络每次广播由两个邻居接收
Etx_mJ = 0.366;
Erx_mJ = 0.310;
activeCommunicationEnergy_J = (Etx_mJ*nBroadcast+Erx_mJ*nReceive)/1000;

m = struct;
m.NormalizedTrackingRMSE = sqrt(mean(z.^2,'all'));
m.NormalizedTrackingRMSE_First6h = sqrt(mean(z(first6,:,:).^2,'all'));
m.NormalizedTrackingRMSE_After6h = sqrt(mean(z(after6,:,:).^2,'all'));
m.TemperatureRMSE_C = sqrt(mean(dT.^2,'all'));
m.HumidityRatioRMSE = sqrt(mean(dW.^2,'all'));
m.CO2RMSE_ppm = sqrt(mean(dC.^2,'all'));
m.GraphDisagreementRMSE = sqrt(mean(graphDisagreement.^2));
m.GraphDisagreementRMSE_First6h = sqrt(mean(graphDisagreement(first6).^2));
m.SameStagePairRMSE = sqrt(mean(sameStageDisagreement.^2));
m.SameStagePairRMSE_First6h = sqrt(mean(sameStageDisagreement(first6).^2));
m.MaxAbsNormalizedError = max(abs(z),[],'all');
m.SafetyViolationSamples = nnz(anyViolation);
m.TemperatureViolationSamples = nnz(violTLow|violTHigh);
m.HumidityViolationSamples = nnz(violWLow|violWHigh);
m.CO2ViolationSamples = nnz(violCLow|violCHigh);
m.SafetyCompliance_pct = 100*(1-nnz(anyViolation)/numel(anyViolation));
m.SaturationAtControlUpdates = ...
    nnz(result.satQ(controlIndex,:))+nnz(result.satW(controlIndex,:)) ...
    +nnz(result.satF(controlIndex,:));
m.DeviceEnergy_kWh = sum(energyPerRoom);
m.ControlExecutions = numel(controlIndex)*nRoom;
m.NetworkTransmissions = nBroadcast;
m.NetworkReceptions = nReceive;
m.ActiveCommunicationEnergy_J = activeCommunicationEnergy_J;

roomRMSET = sqrt(mean(dT.^2,1))';
roomRMSEW = sqrt(mean(dW.^2,1))';
roomRMSEC = sqrt(mean(dC.^2,1))';
roomNormRMSE = squeeze(sqrt(mean(mean(z.^2,1),3)))';
roomViolation = sum(anyViolation,1)';
roomSat = sum(result.satQ(controlIndex,:),1)' ...
    +sum(result.satW(controlIndex,:),1)' ...
    +sum(result.satF(controlIndex,:),1)';

result.Cppm=Cppm; result.dT=dT; result.dW=dW; result.dC=dC;
result.z=z; result.anyViolation=anyViolation;
result.violTLow=violTLow; result.violTHigh=violTHigh;
result.violWLow=violWLow; result.violWHigh=violWHigh;
result.violCLow=violCLow; result.violCHigh=violCHigh;
result.devicePower_kW=devicePower; result.energyPerRoom_kWh=energyPerRoom;
result.graphDisagreement=graphDisagreement;
result.sameStageDisagreement=sameStageDisagreement;
result.trackingTime=trackingTime;
result.metric=m;
result.roomMetric=struct('RMSET',roomRMSET,'RMSEW',roomRMSEW, ...
    'RMSEC',roomRMSEC,'NormRMSE',roomNormRMSE, ...
    'ViolationSamples',roomViolation,'SaturationUpdates',roomSat, ...
    'Energy_kWh',energyPerRoom);
end

function tableOut = build_comparison_table(independent,collaborative)
Condition = {'Independent_gamma_0';'Collaborative_gamma_0_30'};
Gamma = [independent.gamma;collaborative.gamma];
NormalizedTrackingRMSE = [independent.metric.NormalizedTrackingRMSE; ...
    collaborative.metric.NormalizedTrackingRMSE];
NormalizedTrackingRMSE_First6h = ...
    [independent.metric.NormalizedTrackingRMSE_First6h; ...
     collaborative.metric.NormalizedTrackingRMSE_First6h];
NormalizedTrackingRMSE_After6h = ...
    [independent.metric.NormalizedTrackingRMSE_After6h; ...
     collaborative.metric.NormalizedTrackingRMSE_After6h];
TemperatureRMSE_C = [independent.metric.TemperatureRMSE_C; ...
    collaborative.metric.TemperatureRMSE_C];
HumidityRatioRMSE = [independent.metric.HumidityRatioRMSE; ...
    collaborative.metric.HumidityRatioRMSE];
CO2RMSE_ppm = [independent.metric.CO2RMSE_ppm; ...
    collaborative.metric.CO2RMSE_ppm];
GraphDisagreementRMSE = [independent.metric.GraphDisagreementRMSE; ...
    collaborative.metric.GraphDisagreementRMSE];
GraphDisagreementRMSE_First6h = ...
    [independent.metric.GraphDisagreementRMSE_First6h; ...
     collaborative.metric.GraphDisagreementRMSE_First6h];
SameStagePairRMSE = [independent.metric.SameStagePairRMSE; ...
    collaborative.metric.SameStagePairRMSE];
SameStagePairRMSE_First6h = ...
    [independent.metric.SameStagePairRMSE_First6h; ...
     collaborative.metric.SameStagePairRMSE_First6h];
MaxAbsNormalizedError = [independent.metric.MaxAbsNormalizedError; ...
    collaborative.metric.MaxAbsNormalizedError];
SafetyViolationSamples = [independent.metric.SafetyViolationSamples; ...
    collaborative.metric.SafetyViolationSamples];
TemperatureViolationSamples = ...
    [independent.metric.TemperatureViolationSamples; ...
     collaborative.metric.TemperatureViolationSamples];
HumidityViolationSamples = [independent.metric.HumidityViolationSamples; ...
    collaborative.metric.HumidityViolationSamples];
CO2ViolationSamples = [independent.metric.CO2ViolationSamples; ...
    collaborative.metric.CO2ViolationSamples];
SafetyCompliance_pct = [independent.metric.SafetyCompliance_pct; ...
    collaborative.metric.SafetyCompliance_pct];
SaturationAtControlUpdates = ...
    [independent.metric.SaturationAtControlUpdates; ...
     collaborative.metric.SaturationAtControlUpdates];
DeviceEnergy_kWh = [independent.metric.DeviceEnergy_kWh; ...
    collaborative.metric.DeviceEnergy_kWh];
ControlExecutions = [independent.metric.ControlExecutions; ...
    collaborative.metric.ControlExecutions];
NetworkTransmissions = [independent.metric.NetworkTransmissions; ...
    collaborative.metric.NetworkTransmissions];
NetworkReceptions = [independent.metric.NetworkReceptions; ...
    collaborative.metric.NetworkReceptions];
ActiveCommunicationEnergy_J = ...
    [independent.metric.ActiveCommunicationEnergy_J; ...
     collaborative.metric.ActiveCommunicationEnergy_J];
tableOut = table(Condition,Gamma,NormalizedTrackingRMSE, ...
    NormalizedTrackingRMSE_First6h,NormalizedTrackingRMSE_After6h, ...
    TemperatureRMSE_C,HumidityRatioRMSE,CO2RMSE_ppm, ...
    GraphDisagreementRMSE,GraphDisagreementRMSE_First6h, ...
    SameStagePairRMSE,SameStagePairRMSE_First6h, ...
    MaxAbsNormalizedError,SafetyViolationSamples, ...
    TemperatureViolationSamples,HumidityViolationSamples, ...
    CO2ViolationSamples,SafetyCompliance_pct, ...
    SaturationAtControlUpdates,DeviceEnergy_kWh,ControlExecutions, ...
    NetworkTransmissions,NetworkReceptions,ActiveCommunicationEnergy_J);
end

function tableOut = build_effect_table(independent,collaborative)
Metric = {'NormalizedTrackingRMSE';'NormalizedTrackingRMSE_First6h'; ...
    'NormalizedTrackingRMSE_After6h';'GraphDisagreementRMSE'; ...
    'GraphDisagreementRMSE_First6h';'SameStagePairRMSE'; ...
    'SameStagePairRMSE_First6h';'TemperatureRMSE_C'; ...
    'HumidityRatioRMSE';'CO2RMSE_ppm';'DeviceEnergy_kWh'};
Independent = [independent.metric.NormalizedTrackingRMSE; ...
    independent.metric.NormalizedTrackingRMSE_First6h; ...
    independent.metric.NormalizedTrackingRMSE_After6h; ...
    independent.metric.GraphDisagreementRMSE; ...
    independent.metric.GraphDisagreementRMSE_First6h; ...
    independent.metric.SameStagePairRMSE; ...
    independent.metric.SameStagePairRMSE_First6h; ...
    independent.metric.TemperatureRMSE_C; ...
    independent.metric.HumidityRatioRMSE;independent.metric.CO2RMSE_ppm; ...
    independent.metric.DeviceEnergy_kWh];
Collaborative = [collaborative.metric.NormalizedTrackingRMSE; ...
    collaborative.metric.NormalizedTrackingRMSE_First6h; ...
    collaborative.metric.NormalizedTrackingRMSE_After6h; ...
    collaborative.metric.GraphDisagreementRMSE; ...
    collaborative.metric.GraphDisagreementRMSE_First6h; ...
    collaborative.metric.SameStagePairRMSE; ...
    collaborative.metric.SameStagePairRMSE_First6h; ...
    collaborative.metric.TemperatureRMSE_C; ...
    collaborative.metric.HumidityRatioRMSE;collaborative.metric.CO2RMSE_ppm; ...
    collaborative.metric.DeviceEnergy_kWh];
AbsoluteChange_CollabMinusIndependent = Collaborative-Independent;
Improvement_pct_PositiveIsBetter = 100*(Independent-Collaborative)./Independent;
Direction = repmat({'lower_is_better'},numel(Metric),1);
tableOut = table(Metric,Independent,Collaborative, ...
    AbsoluteChange_CollabMinusIndependent, ...
    Improvement_pct_PositiveIsBetter,Direction);
end

function tableOut = build_room_comparison(independent,collaborative, ...
    roomNames,stageNames,roomStage)
Room = [roomNames;roomNames];
Stage = [stageNames(roomStage);stageNames(roomStage)];
Condition = [repmat({'Independent_gamma_0'},6,1); ...
    repmat({'Collaborative_gamma_0_30'},6,1)];
RMSE_T_C = [independent.roomMetric.RMSET;collaborative.roomMetric.RMSET];
RMSE_w_kg_kgda = [independent.roomMetric.RMSEW;collaborative.roomMetric.RMSEW];
RMSE_CO2_ppm = [independent.roomMetric.RMSEC;collaborative.roomMetric.RMSEC];
NormalizedRMSE = [independent.roomMetric.NormRMSE; ...
    collaborative.roomMetric.NormRMSE];
SafetyViolationSamples = [independent.roomMetric.ViolationSamples; ...
    collaborative.roomMetric.ViolationSamples];
SaturationUpdates = [independent.roomMetric.SaturationUpdates; ...
    collaborative.roomMetric.SaturationUpdates];
DeviceEnergy_kWh = [independent.roomMetric.Energy_kWh; ...
    collaborative.roomMetric.Energy_kWh];
tableOut = table(Room,Stage,Condition,RMSE_T_C,RMSE_w_kg_kgda, ...
    RMSE_CO2_ppm,NormalizedRMSE,SafetyViolationSamples, ...
    SaturationUpdates,DeviceEnergy_kWh);
end

function write_case_outputs(resultDir,prefix,result,clock,outdoor, ...
    roomNames,stageNames,roomStage,model)
nStep = clock.nStep;
nRoom = clock.nRoom;
timeLong = repmat(clock.t,nRoom,1);
roomIndexLong = kron((1:nRoom)',ones(nStep,1));
roomLong = roomNames(roomIndexLong);
stageLong = stageNames(roomStage(roomIndexLong));
weatherUpdate = mod(clock.t,3600)==0;

stateTable = table(timeLong,timeLong/3600,roomLong,stageLong, ...
    result.T(:),result.w(:),result.Cppm(:),result.dT(:),result.dW(:), ...
    result.dC(:),reshape(result.z(:,:,1),[],1), ...
    reshape(result.z(:,:,2),[],1),reshape(result.z(:,:,3),[],1), ...
    result.violTLow(:),result.violTHigh(:),result.violWLow(:), ...
    result.violWHigh(:),result.violCLow(:),result.violCHigh(:), ...
    result.anyViolation(:),repmat(result.isControlUpdate,nRoom,1), ...
    result.isBroadcast(:),repmat(weatherUpdate,nRoom,1), ...
    repmat(outdoor.T,nRoom,1),repmat(outdoor.w,nRoom,1), ...
    repmat(outdoor.Cppm,nRoom,1), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','T_C','w_kg_kgda', ...
    'CO2_ppm','dT_C','dw_kg_kgda','dCO2_ppm','zT','zw','zCO2', ...
    'T_Low','T_High','w_Low','w_High','CO2_Low','CO2_High', ...
    'AnySafetyViolation','ControlUpdate','Broadcast','WeatherUpdate', ...
    'T_out_C','w_out_kg_kgda','CO2_out_ppm'});
writetable(stateTable,fullfile(resultDir,[prefix '_state_trajectories.csv']));

[cool,heat,humid,dehumid,fanLow,fanHigh,fanOff,coolDuty,heatDuty, ...
    humidDuty,dehumidDuty,power] = decode_devices( ...
    result.Qcmd,result.mwCmd,result.Fcmd,model.device);
deviceTable = table(timeLong,timeLong/3600,roomLong,stageLong, ...
    result.Qraw(:),result.Qcmd(:),cool(:),heat(:),coolDuty(:),heatDuty(:), ...
    result.mwRaw(:),result.mwCmd(:),humid(:),dehumid(:), ...
    humidDuty(:),dehumidDuty(:),result.Fraw(:),result.Fcmd(:), ...
    result.Fcmd(:)*3600,fanOff(:),fanLow(:),fanHigh(:),power(:), ...
    result.satQ(:),result.satW(:),result.satF(:), ...
    reshape(result.satDetail(:,:,1),[],1), ...
    reshape(result.satDetail(:,:,2),[],1), ...
    reshape(result.satDetail(:,:,3),[],1), ...
    reshape(result.satDetail(:,:,4),[],1), ...
    reshape(result.satDetail(:,:,5),[],1), ...
    reshape(result.satDetail(:,:,6),[],1), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','Q_raw_W','Q_cmd_W', ...
    'Cooling_kW','Heating_kW','CoolingDuty','HeatingDuty','mw_raw_kg_s', ...
    'mw_cmd_kg_s','Humidification_kg_h','Dehumidification_kg_h', ...
    'HumidifierDuty','DehumidifierDuty','F_raw_m3_s','F_cmd_m3_s', ...
    'F_cmd_m3_h','FanOffDuty','FanLowDuty','FanHighDuty','DevicePower_kW', ...
    'ThermalSaturation','MoistureSaturation','VentSaturation', ...
    'CoolingSaturated','HeatingSaturated','DehumidificationSaturated', ...
    'HumidificationSaturated','VentLowSaturated','VentHighSaturated'});
writetable(deviceTable,fullfile(resultDir,[prefix '_device_commands.csv']));

violationTable = stateTable(result.anyViolation(:),:);
writetable(violationTable, ...
    fullfile(resultDir,[prefix '_constraint_violations.csv']));

updateMask = false(nStep,nRoom);
updateMask(result.isControlUpdate,:)=true;
writetable(deviceTable(updateMask(:),:), ...
    fullfile(resultDir,[prefix '_commands_at_control_updates.csv']));
end

function [weather,selection] = select_typical_72h( ...
    date,T,RH,P,CO2ppm,w)
nHour = numel(T);
startIndex = (1:24:(nHour-72+1))';
nCandidate = numel(startIndex);
feature = zeros(nCandidate,4);
for q = 1:nCandidate
    idx = startIndex(q):(startIndex(q)+71);
    feature(q,:) = [mean(T(idx)),mean(w(idx)),std(T(idx),1),std(w(idx),1)];
end
center = [mean(T),mean(w),0,0];
spread = zeros(1,4);
for j = 1:4
    if j>=3
        center(j)=local_quantile(feature(:,j),0.5);
    end
    spread(j)=local_quantile(feature(:,j),0.75) ...
        -local_quantile(feature(:,j),0.25);
end
if any(spread<=0)
    error('典型气象窗口的IQR标准化失败。');
end
score = sum(((feature-center)./spread).^2,2);
[bestScore,bestCandidate] = min(score);
idx = startIndex(bestCandidate):(startIndex(bestCandidate)+71);

weather = struct('Date',date(idx),'T',T(idx),'RH',RH(idx),'P',P(idx), ...
    'CO2ppm',CO2ppm(idx),'w',w(idx));
selection = struct('SelectedStart',date(idx(1)), ...
    'SelectedEnd',date(idx(end)),'CandidateCount',nCandidate, ...
    'Score',bestScore,'MeanT',feature(bestCandidate,1), ...
    'MeanW',feature(bestCandidate,2),'StdT',feature(bestCandidate,3), ...
    'StdW',feature(bestCandidate,4),'MinT',min(T(idx)), ...
    'MaxT',max(T(idx)),'MinRH',min(RH(idx)),'MaxRH',max(RH(idx)), ...
    'StartExcelRow',idx(1)+1);
end

function q = local_quantile(x,p)
x = sort(x(:));
position = 1+(numel(x)-1)*p;
lo = floor(position);
hi = ceil(position);
q = x(lo)+(x(hi)-x(lo))*(position-lo);
end

function w = humidity_ratio_from_t_rh_p(TdegC,RH_pct,P_Pa)
Tk = TdegC+273.15;
lnPws = zeros(size(TdegC));
above = TdegC>=0;
lnPws(above) = -5.8002206e3./Tk(above)+1.3914993 ...
    -4.8640239e-2.*Tk(above)+4.1764768e-5.*Tk(above).^2 ...
    -1.4452093e-8.*Tk(above).^3+6.5459673.*log(Tk(above));
below = ~above;
lnPws(below) = -5.6745359e3./Tk(below)+6.3925247 ...
    -9.677843e-3.*Tk(below)+6.2215701e-7.*Tk(below).^2 ...
    +2.0747825e-9.*Tk(below).^3-9.484024e-13.*Tk(below).^4 ...
    +4.1635019.*log(Tk(below));
Pws = exp(lnPws);
Pv = (RH_pct/100).*Pws;
w = 0.621945*Pv./(P_Pa-Pv);
end

function [cool_kW,heat_kW,humid_kg_h,dehumid_kg_h, ...
    fanLowDuty,fanHighDuty,fanOffDuty,coolDuty,heatDuty, ...
    humidDuty,dehumidDuty,totalPower_kW] = decode_devices(Q,mw,F,device)
cool_kW = max(-Q,0)/1000;
heat_kW = max(Q,0)/1000;
humid_kg_h = max(mw,0)*3600;
dehumid_kg_h = max(-mw,0)*3600;
coolDuty = cool_kW/(device.coolRated_W/1000);
heatDuty = heat_kW/(device.heatRated_W/1000);
humidDuty = humid_kg_h/device.humMax_kg_h;
dehumidDuty = dehumid_kg_h/ ...
    (device.nDehumidifier*device.dehPerUnit_kg_h);

flow = F*3600;
fanLowDuty = zeros(size(F));
fanHighDuty = zeros(size(F));
low = flow<=device.fanLow_m3_h;
fanLowDuty(low)=flow(low)/device.fanLow_m3_h;
high = ~low;
fanHighDuty(high)=(flow(high)-device.fanLow_m3_h) ...
    /(device.fanHigh_m3_h-device.fanLow_m3_h);
fanLowDuty(high)=1-fanHighDuty(high);
fanOffDuty=max(0,1-fanLowDuty-fanHighDuty);

totalPower_kW = device.coolPowerRated_kW*coolDuty ...
    +device.heatPowerRated_kW*heatDuty ...
    +device.humPowerRated_kW*humidDuty ...
    +(device.nDehumidifier*device.dehPowerPerUnit_kW)*dehumidDuty ...
    +device.fanLowPower_kW*fanLowDuty ...
    +device.fanHighPower_kW*fanHighDuty;
end

function write_provenance(filePath,weatherFile,selection,dt,controlPeriod, ...
    weatherPeriod,tEnd,gamma0,gamma1,device,limits)
Parameter = {'weatherFile';'selectedWeatherStart';'weatherSelectionRule'; ...
    'dt';'controlPeriod';'weatherPeriod';'tEnd';'gammaIndependent'; ...
    'gammaCollaborative';'topology';'initialRoom2HumidityRatio'; ...
    'coolingLimit';'heatingLimit';'humidificationLimit'; ...
    'dehumidificationUnits';'dehumidificationLimit';'ventilationLimit'; ...
    'dynamicWeatherFeedforward';'trackingScale';'controlGainH'; ...
    'activeTxEnergy';'activeRxEnergy'};
Value = {weatherFile;char(selection.SelectedStart); ...
    ['72 h daily windows; means target annual means; std targets candidate ', ...
     'medians; candidate IQR scaling; equal weights']; ...
    sprintf('%g',dt);sprintf('%g',controlPeriod);sprintf('%g',weatherPeriod); ...
    sprintf('%g',tEnd);sprintf('%.2f',gamma0);sprintf('%.2f',gamma1); ...
    'undirected ring 1-2-3-4-5-6-1';'0.0165089'; ...
    sprintf('%g',limits.QMin);sprintf('%g',limits.QMax); ...
    sprintf('%g',limits.mwMax*3600);sprintf('%d',device.nDehumidifier); ...
    sprintf('%g',-limits.mwMin*3600);sprintf('%g',limits.FMax*3600); ...
    'step1 equilibrium equations recomputed from hourly outdoor state'; ...
    'target-to-nearest-safety-bound';'same as experiment1';'0.366';'0.310'};
Unit = {'path';'datetime';'-';'s';'s';'s';'s';'-';'-';'-'; ...
    'kg/kgda';'W';'W';'kg/h';'count';'kg/h';'m3/h';'-';'-';'-'; ...
    'mJ/transmission';'mJ/reception'};
Category = {'用户提供文件';'派生计算值';'算法设计值';'用户指定值'; ...
    '用户指定值';'用户指定值';'用户指定值';'算法设计值'; ...
    '算法设计值';'算法设计值';'用户指定值';'厂商原始规格'; ...
    '厂商原始规格';'厂商原始规格';'算法/设备配置值'; ...
    '派生计算值';'厂商原始规格';'派生计算值';'派生计算值'; ...
    '算法设计值';'文献/厂商参数派生值';'文献/厂商参数派生值'};
SourceNote = {'用户提供的南宁2023逐时重建气象表'; ...
    '按固定典型窗口规则自动得到，非人工挑选'; ...
    '消融实验统一的典型动态气象选择规则'; ...
    '与实验一一致';'与实验一一致';'与用户要求一致'; ...
    '与用户要求一致';'删除协同项';'与实验一一致'; ...
    '与实验一一致';'明确排除作废值0.0171547'; ...
    '继承实验一';'继承实验一';'继承实验一';'继承实验一'; ...
    '3台乘单台0.8375 kg/h';'继承实验一'; ...
    '由step1平衡方程随相同逐时气象派生，两组一致'; ...
    '由阶段安全边界派生';'与实验一相同H_g'; ...
    '代表性CC2530参数，仅估算主动通信'; ...
    '代表性CC2530参数，仅估算主动通信'};
tableOut = table(Parameter,Value,Unit,Category,SourceNote);
writetable(tableOut,filePath);
end

function make_figures(resultDir,clock,weather,independent,collaborative,roomNames)
timeHour = clock.t/3600;
fontName = choose_chinese_font();
figureVisibility = 'on'; % 改为'off'时只保存图片，不弹出图窗

fig1=figure('Visible',figureVisibility,'Color','w','Position',[100 80 1120 780]);
subplot(4,1,1); stairs((0:71),weather.T,'LineWidth',1.2); grid on;
ylabel('室外温度/℃'); title('选取的典型72 h室外气象');
subplot(4,1,2); stairs((0:71),weather.RH,'LineWidth',1.2); grid on;
ylabel('室外相对湿度/%');
subplot(4,1,3); stairs((0:71),weather.w,'LineWidth',1.2); grid on;
ylabel('室外含湿量/(kg/kg干空气)');
subplot(4,1,4); stairs((0:71),weather.P/1000,'LineWidth',1.2); grid on;
ylabel('室外压力/kPa'); xlabel('时间/h');
set(findall(fig1,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig1,fullfile(resultDir,'selected_weather_72h.png'), ...
    'Resolution',150);

fig2=figure('Visible',figureVisibility,'Color','w','Position',[120 100 1120 760]);
subplot(3,1,1);
plot(timeHour,independent.trackingTime,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.trackingTime,'LineWidth',1.2);
grid on; ylabel('综合归一化跟踪误差');
title('独立控制与协同控制的跟踪及一致性对比');
legend({'独立控制','协同控制'});
subplot(3,1,2);
plot(timeHour,independent.graphDisagreement,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.graphDisagreement,'LineWidth',1.2);
grid on; ylabel('图一致性误差');
subplot(3,1,3);
plot(timeHour,independent.sameStageDisagreement,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.sameStageDisagreement,'LineWidth',1.2);
grid on; ylabel('同阶段仓对误差'); xlabel('时间/h');
set(findall(fig2,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig2,fullfile(resultDir,'tracking_and_disagreement.png'), ...
    'Resolution',150);

fig2b=figure('Visible',figureVisibility,'Color','w','Position',[130 110 1120 760]);
subplot(3,1,1);
plot(timeHour,independent.trackingTime,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.trackingTime,'LineWidth',1.2);
grid on; xlim([0 6]); ylabel('综合归一化跟踪误差');
title('前6 h跟踪及一致性对比');
legend({'独立控制','协同控制'});
subplot(3,1,2);
plot(timeHour,independent.graphDisagreement,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.graphDisagreement,'LineWidth',1.2);
grid on; xlim([0 6]); ylabel('图一致性误差');
subplot(3,1,3);
plot(timeHour,independent.sameStageDisagreement,'LineWidth',1.2); hold on;
plot(timeHour,collaborative.sameStageDisagreement,'LineWidth',1.2);
grid on; xlim([0 6]); ylabel('同阶段仓对误差'); xlabel('时间/h');
set(findall(fig2b,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig2b, ...
    fullfile(resultDir,'tracking_and_disagreement_first6h.png'), ...
    'Resolution',150);

fig3=figure('Visible',figureVisibility,'Color','w','Position',[140 120 1160 760]);
pairs=[1 2;4 5];
for p=1:2
    subplot(2,1,p);
    localDiff=sqrt(mean((squeeze(independent.z(:,pairs(p,1),:)) ...
        -squeeze(independent.z(:,pairs(p,2),:))).^2,2));
    collabDiff=sqrt(mean((squeeze(collaborative.z(:,pairs(p,1),:)) ...
        -squeeze(collaborative.z(:,pairs(p,2),:))).^2,2));
    plot(timeHour,localDiff,'LineWidth',1.2); hold on;
    plot(timeHour,collabDiff,'LineWidth',1.2); grid on;
    ylabel('归一化偏差差异');
    title(sprintf('%s—%s同阶段仓对',roomNames{pairs(p,1)}, ...
        roomNames{pairs(p,2)}));
    legend({'独立控制','协同控制'});
end
xlabel('时间/h');
set(findall(fig3,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig3,fullfile(resultDir,'same_stage_pair_comparison.png'), ...
    'Resolution',150);

fig4=figure('Visible',figureVisibility,'Color','w','Position',[160 140 1140 460]);
values=[independent.metric.NormalizedTrackingRMSE, ...
        collaborative.metric.NormalizedTrackingRMSE; ...
        independent.metric.GraphDisagreementRMSE, ...
        collaborative.metric.GraphDisagreementRMSE; ...
        independent.metric.SameStagePairRMSE, ...
        collaborative.metric.SameStagePairRMSE];
bar(values); grid on;
set(gca,'XTickLabel',{'综合跟踪','图一致性','同阶段仓对'});
ylabel('归一化RMSE'); legend({'独立控制','协同控制'});
title('协同作用消融实验核心指标（越小越好）');
set(findall(fig4,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig4,fullfile(resultDir,'ablation_metric_bars.png'), ...
    'Resolution',150);

fig5=figure('Visible',figureVisibility,'Color','w','Position',[180 160 1120 720]);
subplot(3,1,1);
plot(timeHour,independent.T,'LineWidth',0.9); grid on;
ylabel('温度/℃'); title('独立控制：六仓状态轨迹');
legend(roomNames,'Location','eastoutside');
subplot(3,1,2); plot(timeHour,independent.w,'LineWidth',0.9); grid on;
ylabel('含湿量/(kg/kg干空气)');
subplot(3,1,3); plot(timeHour,independent.Cppm,'LineWidth',0.9); grid on;
ylabel('CO₂浓度/ppm'); xlabel('时间/h');
set(findall(fig5,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig5, ...
    fullfile(resultDir,'independent_state_trajectories.png'), ...
    'Resolution',150);

fig6=figure('Visible',figureVisibility,'Color','w','Position',[200 180 1120 720]);
subplot(3,1,1);
plot(timeHour,collaborative.T,'LineWidth',0.9); grid on;
ylabel('温度/℃'); title('协同控制：六仓状态轨迹');
legend(roomNames,'Location','eastoutside');
subplot(3,1,2); plot(timeHour,collaborative.w,'LineWidth',0.9); grid on;
ylabel('含湿量/(kg/kg干空气)');
subplot(3,1,3); plot(timeHour,collaborative.Cppm,'LineWidth',0.9); grid on;
ylabel('CO₂浓度/ppm'); xlabel('时间/h');
set(findall(fig6,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig6,fullfile(resultDir,'collaborative_state_trajectories.png'), ...
    'Resolution',150);
drawnow;
end

function fontName = choose_chinese_font
% 优先选用常见中文字体；若系统没有，则保留MATLAB默认字体。
fontName = get(groot,'defaultAxesFontName');
availableFonts = listfonts;
preferredFonts = {'PingFang SC','Microsoft YaHei','SimHei', ...
    'Noto Sans CJK SC','Heiti SC','Songti SC','Arial Unicode MS'};
for k = 1:numel(preferredFonts)
    if any(strcmpi(availableFonts,preferredFonts{k}))
        fontName = preferredFonts{k};
        return;
    end
end
end

function c = ppm_to_mass(Cppm,TdegC,pressurePa,par)
c = Cppm.*1e-6.*pressurePa.*par.MCO2 ...
    ./(par.R.*(TdegC+273.15));
end

function Cppm = mass_to_ppm(c,TdegC,pressurePa,par)
Cppm = c.*par.R.*(TdegC+273.15) ...
    ./(pressurePa.*par.MCO2).*1e6;
end
