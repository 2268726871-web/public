function experiment2_1_collaboration_gain_sensitivity_4units()
%EXPERIMENT2_1_COLLABORATION_GAIN_SENSITIVITY
% 实验2.1：六仓协同增益gamma敏感性分析。
%
% 目的：
%   在不改变物理模型、控制器H_g、设备能力、动态气象和初始状态的
%   前提下，只扫描协同增益gamma，确定后续周期通信、固定事件触发和
%   自适应事件触发共同使用的协同增益。
%
% 运行：
%   experiment2_1_collaboration_gain_sensitivity
%
% 结果目录：
%   本文件所在目录/results/experiment2_1_gamma_sensitivity/
%
% 数值性质：
%   gamma候选网格、1%跟踪容差和90%协调改善比例均为算法设计值，
%   不是文献原始数据。程序会写入parameter_provenance.csv。

clc;
close all;

%% 1. 仿真时钟和增益候选
dt = 60;                         % s，用户指定
controlPeriod = 300;             % s，用户指定
weatherPeriod = 3600;            % s，用户指定
tEnd = 72*3600;                  % s，用户指定
t = (0:dt:tEnd)';
nStep = numel(t);
nRoom = 6;
nControlStep = controlPeriod/dt;

% [A] 协同增益候选网格。
% 原扫描在gamma=0.30处达到上边界，因此增加0.35和0.40，
% 用于确认0.30是否由1%跟踪性能约束形成合理选择，而非区间截断结果。
gammaGrid = (0:0.05:0.40)';

% [A] 最终增益选择规则中的设计阈值。
trackingTolerance_pct = 1.0;
nearMaximumCoordinationFraction = 0.90;

if mod(controlPeriod,dt)~=0 || mod(weatherPeriod,dt)~=0
    error('控制周期和气象周期必须是积分步长的整数倍。');
end

weatherFile = ...
    '/Users/bystander/Desktop/蘑菇仓/南宁2023逐时室外气象与CO2背景数据_重建版.xlsx';
if ~isfile(weatherFile)
    error(['未找到气象文件：\n%s\n请修改weatherFile。'],weatherFile);
end

%% 2. 物理模型参数（继承step1、实验一和实验二）
par.rhoAir = 1.2;
par.cpAir = 1006;
par.UA = 31.5;
par.CT = 6.85e6;
par.V = 104;
par.R = 8.314462618;
par.MCO2 = 44.01e-3;
par.PIndoorRef = 101325;

%% 3. 读取气象并使用与实验二相同的典型72 h
raw = readmatrix(weatherFile,'Sheet','第四部分模型输入', ...
    'Range','A2:G8761');
if size(raw,1)~=8760 || size(raw,2)<6
    error('气象表尺寸应为8760行且至少6列。');
end
if any(~isfinite(raw(:,1:6)),'all')
    error('气象表A:F中存在缺失或非数值数据。');
end

allDate = datetime(raw(:,1),'ConvertFrom','excel');
allTout = raw(:,2);
allRH = raw(:,3);
allP = raw(:,5);
allCO2 = raw(:,6);
allWout = humidity_ratio_from_t_rh_p(allTout,allRH,allP);

[weather,selection] = select_typical_72h( ...
    allDate,allTout,allRH,allP,allCO2,allWout);

weatherHourIndex = min(floor(t/weatherPeriod)+1,72);
outdoor.T = weather.T(weatherHourIndex);
outdoor.w = weather.w(weatherHourIndex);
outdoor.P = weather.P(weatherHourIndex);
outdoor.Cppm = weather.CO2ppm(weatherHourIndex);
outdoor.c = ppm_to_mass( ...
    outdoor.Cppm,outdoor.T,outdoor.P,par);

%% 4. 阶段目标、生物扰动和六仓配置
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
% 2号仓使用有效值0.0165089，不使用作废值0.0171547。
initial.w = [0.0166165,0.0165089,0.0175666, ...
             0.0098636,0.0094054,0.0099782];
initial.Cppm = [13500,13000,13750,960,930,940];

%% 5. 设备能力、拓扑和反馈增益
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

model = struct('par',par,'T0',T0,'w0',w0,'Cppm0',Cppm0, ...
    'Mda',Mda,'Qbio',Qbio,'mbioW',mbioW,'mbioC',mbioC, ...
    'TMin',TMin,'TMax',TMax,'wMin',wMin,'wMax',wMax, ...
    'CppmMin',CppmMin,'CppmMax',CppmMax,'scaleT',scaleT, ...
    'scaleW',scaleW,'scaleCppm',scaleCppm,'roomStage',roomStage, ...
    'H',H,'L',L,'A',A,'limits',limits,'device',device);
clock = struct('dt',dt,'controlPeriod',controlPeriod, ...
    'nControlStep',nControlStep,'t',t,'nStep',nStep,'nRoom',nRoom);

%% 6. 逐个扫描gamma
nGamma = numel(gammaGrid);
caseResults = cell(nGamma,1);

NormalizedTrackingRMSE = zeros(nGamma,1);
NormalizedTrackingRMSE_First6h = zeros(nGamma,1);
NormalizedTrackingRMSE_After6h = zeros(nGamma,1);
TemperatureRMSE_C = zeros(nGamma,1);
HumidityRatioRMSE = zeros(nGamma,1);
CO2RMSE_ppm = zeros(nGamma,1);
GraphDisagreementRMSE = zeros(nGamma,1);
GraphDisagreementRMSE_First6h = zeros(nGamma,1);
SameStagePairRMSE = zeros(nGamma,1);
SameStagePairRMSE_First6h = zeros(nGamma,1);
SafetyViolationSamples = zeros(nGamma,1);
TemperatureViolationSamples = zeros(nGamma,1);
HumidityViolationSamples = zeros(nGamma,1);
CO2ViolationSamples = zeros(nGamma,1);
SaturationAtControlUpdates = zeros(nGamma,1);
DeviceEnergy_kWh = zeros(nGamma,1);
MaxAbsNormalizedError = zeros(nGamma,1);

roomNormRMSE = zeros(nGamma,nRoom);
roomRMSET = zeros(nGamma,nRoom);
roomRMSEW = zeros(nGamma,nRoom);
roomRMSEC = zeros(nGamma,nRoom);
roomViolation = zeros(nGamma,nRoom);
roomSaturation = zeros(nGamma,nRoom);

fprintf('\n开始协同增益敏感性扫描：\n');
for q = 1:nGamma
    gamma = gammaGrid(q);
    fprintf('  [%d/%d] gamma = %.2f\n',q,nGamma,gamma);

    result = simulate_gamma_case(gamma,initial,model,clock,outdoor);
    result = analyze_gamma_case(result,model,clock);
    caseResults{q} = result;

    m = result.metric;
    NormalizedTrackingRMSE(q) = m.NormalizedTrackingRMSE;
    NormalizedTrackingRMSE_First6h(q) = ...
        m.NormalizedTrackingRMSE_First6h;
    NormalizedTrackingRMSE_After6h(q) = ...
        m.NormalizedTrackingRMSE_After6h;
    TemperatureRMSE_C(q) = m.TemperatureRMSE_C;
    HumidityRatioRMSE(q) = m.HumidityRatioRMSE;
    CO2RMSE_ppm(q) = m.CO2RMSE_ppm;
    GraphDisagreementRMSE(q) = m.GraphDisagreementRMSE;
    GraphDisagreementRMSE_First6h(q) = ...
        m.GraphDisagreementRMSE_First6h;
    SameStagePairRMSE(q) = m.SameStagePairRMSE;
    SameStagePairRMSE_First6h(q) = ...
        m.SameStagePairRMSE_First6h;
    SafetyViolationSamples(q) = m.SafetyViolationSamples;
    TemperatureViolationSamples(q) = m.TemperatureViolationSamples;
    HumidityViolationSamples(q) = m.HumidityViolationSamples;
    CO2ViolationSamples(q) = m.CO2ViolationSamples;
    SaturationAtControlUpdates(q) = m.SaturationAtControlUpdates;
    DeviceEnergy_kWh(q) = m.DeviceEnergy_kWh;
    MaxAbsNormalizedError(q) = m.MaxAbsNormalizedError;

    roomNormRMSE(q,:) = result.roomMetric.NormRMSE';
    roomRMSET(q,:) = result.roomMetric.RMSET';
    roomRMSEW(q,:) = result.roomMetric.RMSEW';
    roomRMSEC(q,:) = result.roomMetric.RMSEC';
    roomViolation(q,:) = result.roomMetric.ViolationSamples';
    roomSaturation(q,:) = result.roomMetric.SaturationUpdates';
end

%% 7. 固定、可复现的gamma选择规则
baseline = 1; % gamma=0
TrackingChange_pct = 100*(NormalizedTrackingRMSE ...
    /NormalizedTrackingRMSE(baseline)-1);
GraphImprovement_pct = 100*(1-GraphDisagreementRMSE ...
    /GraphDisagreementRMSE(baseline));
SameStageImprovement_pct = 100*(1-SameStagePairRMSE ...
    /SameStagePairRMSE(baseline));

CoordinationIndex = 0.5*( ...
    GraphDisagreementRMSE/GraphDisagreementRMSE(baseline) ...
    +SameStagePairRMSE/SameStagePairRMSE(baseline));
CoordinationImprovement_pct = 100*(1-CoordinationIndex);

TrackingFeasible = TrackingChange_pct<=trackingTolerance_pct;
SafetyFeasible = SafetyViolationSamples==0;
SaturationFeasible = ...
    SaturationAtControlUpdates==0;
PositiveCollaborativeGain = gammaGrid>0;
Feasible = TrackingFeasible & SafetyFeasible & ...
    SaturationFeasible & PositiveCollaborativeGain;

if ~any(Feasible)
    error(['没有正gamma同时满足跟踪、安全和饱和约束。', ...
        '请检查gamma_sensitivity_summary.csv，不要强行选择。']);
end

maximumCoordinationImprovement = max(CoordinationImprovement_pct(Feasible));
requiredCoordinationImprovement = ...
    nearMaximumCoordinationFraction*maximumCoordinationImprovement;

candidateIndex = find(Feasible & ...
    CoordinationImprovement_pct>=requiredCoordinationImprovement);
selectedIndex = candidateIndex(1); % gammaGrid已升序，取最小正gamma
selectedGamma = gammaGrid(selectedIndex);
selectedAtUpperBoundary = selectedIndex==nGamma;
Selected = false(nGamma,1);
Selected(selectedIndex) = true;

%% 8. 汇总表和选择结果
summaryTable = table(gammaGrid,NormalizedTrackingRMSE, ...
    NormalizedTrackingRMSE_First6h,NormalizedTrackingRMSE_After6h, ...
    TrackingChange_pct,TemperatureRMSE_C,HumidityRatioRMSE, ...
    CO2RMSE_ppm,GraphDisagreementRMSE, ...
    GraphDisagreementRMSE_First6h,GraphImprovement_pct, ...
    SameStagePairRMSE,SameStagePairRMSE_First6h, ...
    SameStageImprovement_pct,CoordinationIndex, ...
    CoordinationImprovement_pct,SafetyViolationSamples, ...
    TemperatureViolationSamples,HumidityViolationSamples, ...
    CO2ViolationSamples,SaturationAtControlUpdates, ...
    DeviceEnergy_kWh,MaxAbsNormalizedError,TrackingFeasible, ...
    SafetyFeasible,SaturationFeasible,Feasible,Selected, ...
    'VariableNames',{'Gamma','NormalizedTrackingRMSE', ...
    'NormalizedTrackingRMSE_First6h','NormalizedTrackingRMSE_After6h', ...
    'TrackingChange_pct_vs_gamma0','TemperatureRMSE_C', ...
    'HumidityRatioRMSE','CO2RMSE_ppm','GraphDisagreementRMSE', ...
    'GraphDisagreementRMSE_First6h','GraphImprovement_pct_vs_gamma0', ...
    'SameStagePairRMSE','SameStagePairRMSE_First6h', ...
    'SameStageImprovement_pct_vs_gamma0','CoordinationIndex', ...
    'CoordinationImprovement_pct_vs_gamma0','SafetyViolationSamples', ...
    'TemperatureViolationSamples','HumidityViolationSamples', ...
    'CO2ViolationSamples','SaturationAtControlUpdates', ...
    'DeviceEnergy_kWh','MaxAbsNormalizedError','TrackingFeasible', ...
    'SafetyFeasible','SaturationFeasible','Feasible','Selected'});

selectionTable = table(selectedGamma,selectedAtUpperBoundary, ...
    trackingTolerance_pct, ...
    nearMaximumCoordinationFraction,maximumCoordinationImprovement, ...
    requiredCoordinationImprovement, ...
    NormalizedTrackingRMSE(selectedIndex), ...
    TrackingChange_pct(selectedIndex), ...
    GraphDisagreementRMSE(selectedIndex), ...
    GraphImprovement_pct(selectedIndex), ...
    SameStagePairRMSE(selectedIndex), ...
    SameStageImprovement_pct(selectedIndex), ...
    SafetyViolationSamples(selectedIndex), ...
    SaturationAtControlUpdates(selectedIndex), ...
    DeviceEnergy_kWh(selectedIndex), ...
    'VariableNames',{'SelectedGamma','SelectedAtUpperBoundary', ...
    'TrackingTolerance_pct', ...
    'NearMaximumCoordinationFraction', ...
    'MaximumFeasibleCoordinationImprovement_pct', ...
    'RequiredCoordinationImprovement_pct', ...
    'SelectedTrackingRMSE','SelectedTrackingChange_pct_vs_gamma0', ...
    'SelectedGraphDisagreementRMSE', ...
    'SelectedGraphImprovement_pct_vs_gamma0', ...
    'SelectedSameStagePairRMSE', ...
    'SelectedSameStageImprovement_pct_vs_gamma0', ...
    'SelectedSafetyViolationSamples', ...
    'SelectedSaturationAtControlUpdates', ...
    'SelectedDeviceEnergy_kWh'});

GammaLong = repelem(gammaGrid,nRoom);
RoomLong = repmat(roomNames,nGamma,1);
StageLong = repmat(stageNames(roomStage),nGamma,1);
roomTable = table(GammaLong,RoomLong,StageLong, ...
    reshape(roomNormRMSE',[],1),reshape(roomRMSET',[],1), ...
    reshape(roomRMSEW',[],1),reshape(roomRMSEC',[],1), ...
    reshape(roomViolation',[],1),reshape(roomSaturation',[],1), ...
    'VariableNames',{'Gamma','Room','Stage','NormalizedRMSE', ...
    'RMSE_T_C','RMSE_w_kg_kgda','RMSE_CO2_ppm', ...
    'SafetyViolationSamples','SaturationUpdates'});

%% 9. 保存CSV、MAT和中文图形
thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(thisDir,'results', ...
    'experiment2_1_gamma_sensitivity_4units');
if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

writetable(summaryTable, ...
    fullfile(resultDir,'gamma_sensitivity_summary.csv'));
writetable(selectionTable, ...
    fullfile(resultDir,'gamma_selection_result.csv'));
writetable(roomTable, ...
    fullfile(resultDir,'gamma_room_metrics.csv'));

weatherTable = table((0:71)',weather.Date,weather.T,weather.RH, ...
    weather.P,weather.w,weather.CO2ppm, ...
    'VariableNames',{'HourIndex','DateTime_Beijing','T_out_C', ...
    'RH_out_pct','P_out_Pa','w_out_kg_kgda','CO2_background_ppm'});
writetable(weatherTable, ...
    fullfile(resultDir,'selected_weather_72h.csv'));

write_provenance(fullfile(resultDir,'parameter_provenance.csv'), ...
    weatherFile,selection,gammaGrid,trackingTolerance_pct, ...
    nearMaximumCoordinationFraction,dt,controlPeriod,weatherPeriod,tEnd);

save(fullfile(resultDir,'gamma_sensitivity.mat'), ...
    'gammaGrid','caseResults','summaryTable','selectionTable', ...
    'roomTable','selectedGamma','selectedIndex', ...
    'selectedAtUpperBoundary','weather','selection', ...
    'model','clock','outdoor');

try
    make_sensitivity_figures(resultDir,gammaGrid,summaryTable, ...
        selectedGamma,trackingTolerance_pct,selectedAtUpperBoundary);
catch ME
    warning('数值结果已保存，但当前图形环境导出PNG失败：%s',ME.message);
end

%% 10. 命令行结果
fprintf('\n============================================================\n');
fprintf('实验2.1：协同增益敏感性分析完成\n');
fprintf('气象：%s 至 %s\n', ...
    char(selection.SelectedStart),char(selection.SelectedEnd));
fprintf('候选gamma：');
fprintf(' %.2f',gammaGrid);
fprintf('\n');
fprintf('选择规则：跟踪增幅<=%.2f%%，安全/饱和不劣于gamma=0，\n', ...
    trackingTolerance_pct);
fprintf('          达到可行候选最大协调改善的至少%.0f%%，取最小gamma。\n', ...
    100*nearMaximumCoordinationFraction);
fprintf('当前候选区间内的推荐值：gamma = %.2f\n',selectedGamma);
if selectedAtUpperBoundary
    fprintf(['注意：推荐值位于扫描区间上限，只能视为区间内结果；', ...
        '正式固定参数前应扩展上限复核。\n']);
end
fprintf('结果目录：%s\n',resultDir);
fprintf('============================================================\n');
disp(summaryTable);
disp(selectionTable);
end

function result = simulate_gamma_case(gamma,initial,model,clock,outdoor)
% 相同仿真函数扫描全部gamma，保证除协同增益外没有隐藏差异。
nStep = clock.nStep;
nRoom = clock.nRoom;
dt = clock.dt;

T = zeros(nStep,nRoom);
w = zeros(nStep,nRoom);
c = zeros(nStep,nRoom);
T(1,:) = initial.T;
w(1,:) = initial.w;
for i = 1:nRoom
    c(1,i) = ppm_to_mass(initial.Cppm(i),T(1,i), ...
        model.par.PIndoorRef,model.par);
end

Qcmd = zeros(nStep,nRoom);
mwCmd = zeros(nStep,nRoom);
Fcmd = zeros(nStep,nRoom);
satQ = false(nStep,nRoom);
satW = false(nStep,nRoom);
satF = false(nStep,nRoom);
isControlUpdate = false(nStep,1);

uHold = zeros(nRoom,3);
satHold = false(nRoom,6);

for k = 1:(nStep-1)
    out = struct('T',outdoor.T(k),'w',outdoor.w(k), ...
        'P',outdoor.P(k),'Cppm',outdoor.Cppm(k),'c',outdoor.c(k));

    if mod(k-1,clock.nControlStep)==0
        % 控制律：
        % z_i为本仓相对所属阶段目标的归一化误差；
        % eta_i=z_i+gamma*sum_j(z_i-z_j)；
        % u_i=u_ff,g(Tout,wout,Pout,Cout)+H_g*eta_i；
        % 最后对热量、湿量和通风指令执行设备限幅。
        [~,uHold,satHold] = controller_step( ...
            T(k,:),w(k,:),c(k,:),model,gamma,out);
        isControlUpdate(k) = true;
    end

    Qcmd(k,:) = uHold(:,1)';
    mwCmd(k,:) = uHold(:,2)';
    Fcmd(k,:) = uHold(:,3)';
    satQ(k,:) = (satHold(:,1)|satHold(:,2))';
    satW(k,:) = (satHold(:,3)|satHold(:,4))';
    satF(k,:) = (satHold(:,5)|satHold(:,6))';

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

out = struct('T',outdoor.T(end),'w',outdoor.w(end), ...
    'P',outdoor.P(end),'Cppm',outdoor.Cppm(end),'c',outdoor.c(end));
[~,uHold,satHold] = controller_step( ...
    T(end,:),w(end,:),c(end,:),model,gamma,out);
isControlUpdate(end) = true;
Qcmd(end,:) = uHold(:,1)';
mwCmd(end,:) = uHold(:,2)';
Fcmd(end,:) = uHold(:,3)';
satQ(end,:) = (satHold(:,1)|satHold(:,2))';
satW(end,:) = (satHold(:,3)|satHold(:,4))';
satF(end,:) = (satHold(:,5)|satHold(:,6))';

result = struct('gamma',gamma,'T',T,'w',w,'c',c, ...
    'Qcmd',Qcmd,'mwCmd',mwCmd,'Fcmd',Fcmd, ...
    'satQ',satQ,'satW',satW,'satF',satF, ...
    'isControlUpdate',isControlUpdate);
end

function [uRaw,uApplied,sat] = controller_step( ...
    Trow,wrow,crow,model,gamma,out)
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

function result = analyze_gamma_case(result,model,clock)
nStep = clock.nStep;
nRoom = clock.nRoom;
Cppm = mass_to_ppm(result.c,result.T, ...
    model.par.PIndoorRef,model.par);

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
    targetT(i)=model.T0(g);
    targetW(i)=model.w0(g);
    targetC(i)=model.Cppm0(g);
    safeTMin(i)=model.TMin(g);
    safeTMax(i)=model.TMax(g);
    safeWMin(i)=model.wMin(g);
    safeWMax(i)=model.wMax(g);
    safeCMin(i)=model.CppmMin(g);
    safeCMax(i)=model.CppmMax(g);
    normT(i)=model.scaleT(g);
    normW(i)=model.scaleW(g);
    normC(i)=model.scaleCppm(g);
end

dT = result.T-targetT;
dW = result.w-targetW;
dC = Cppm-targetC;
zT = dT./normT;
zW = dW./normW;
zC = dC./normC;
z = cat(3,zT,zW,zC);

violT = result.T<safeTMin | result.T>safeTMax;
violW = result.w<safeWMin | result.w>safeWMax;
violC = Cppm<safeCMin | Cppm>safeCMax;
anyViolation = violT|violW|violC;

[~,~,~,~,~,~,~,~,~,~,~,devicePower] = decode_devices( ...
    result.Qcmd,result.mwCmd,result.Fcmd,model.device);
energyPerRoom = sum(devicePower(1:end-1,:),1)'*clock.dt/3600;

edge = [1 2;2 3;3 4;4 5;5 6;6 1];
graphDisagreement = zeros(nStep,1);
for k = 1:nStep
    accumulator = [];
    for e = 1:size(edge,1)
        dz = squeeze(z(k,edge(e,1),:)-z(k,edge(e,2),:));
        accumulator = [accumulator;dz(:)]; %#ok<AGROW>
    end
    graphDisagreement(k) = sqrt(mean(accumulator.^2));
end

sameStagePair = [1 2;4 5];
sameStageDisagreement = zeros(nStep,1);
for k = 1:nStep
    accumulator = [];
    for p = 1:size(sameStagePair,1)
        dz = squeeze(z(k,sameStagePair(p,1),:) ...
            -z(k,sameStagePair(p,2),:));
        accumulator = [accumulator;dz(:)]; %#ok<AGROW>
    end
    sameStageDisagreement(k) = sqrt(mean(accumulator.^2));
end

first6 = clock.t<=6*3600;
after6 = clock.t>=6*3600;
controlIndex = find(result.isControlUpdate);

m = struct;
m.NormalizedTrackingRMSE = sqrt(mean(z.^2,'all'));
m.NormalizedTrackingRMSE_First6h = ...
    sqrt(mean(z(first6,:,:).^2,'all'));
m.NormalizedTrackingRMSE_After6h = ...
    sqrt(mean(z(after6,:,:).^2,'all'));
m.TemperatureRMSE_C = sqrt(mean(dT.^2,'all'));
m.HumidityRatioRMSE = sqrt(mean(dW.^2,'all'));
m.CO2RMSE_ppm = sqrt(mean(dC.^2,'all'));
m.GraphDisagreementRMSE = sqrt(mean(graphDisagreement.^2));
m.GraphDisagreementRMSE_First6h = ...
    sqrt(mean(graphDisagreement(first6).^2));
m.SameStagePairRMSE = sqrt(mean(sameStageDisagreement.^2));
m.SameStagePairRMSE_First6h = ...
    sqrt(mean(sameStageDisagreement(first6).^2));
m.SafetyViolationSamples = nnz(anyViolation);
m.TemperatureViolationSamples = nnz(violT);
m.HumidityViolationSamples = nnz(violW);
m.CO2ViolationSamples = nnz(violC);
m.SaturationAtControlUpdates = ...
    nnz(result.satQ(controlIndex,:)) ...
    +nnz(result.satW(controlIndex,:)) ...
    +nnz(result.satF(controlIndex,:));
m.DeviceEnergy_kWh = sum(energyPerRoom);
m.MaxAbsNormalizedError = max(abs(z),[],'all');

roomNormRMSE = squeeze(sqrt(mean(mean(z.^2,1),3)))';
roomRMSET = sqrt(mean(dT.^2,1))';
roomRMSEW = sqrt(mean(dW.^2,1))';
roomRMSEC = sqrt(mean(dC.^2,1))';
roomViolation = sum(anyViolation,1)';
roomSaturation = sum(result.satQ(controlIndex,:),1)' ...
    +sum(result.satW(controlIndex,:),1)' ...
    +sum(result.satF(controlIndex,:),1)';

result.Cppm = Cppm;
result.z = z;
result.graphDisagreement = graphDisagreement;
result.sameStageDisagreement = sameStageDisagreement;
result.metric = m;
result.roomMetric = struct('NormRMSE',roomNormRMSE, ...
    'RMSET',roomRMSET,'RMSEW',roomRMSEW,'RMSEC',roomRMSEC, ...
    'ViolationSamples',roomViolation, ...
    'SaturationUpdates',roomSaturation);
end

function [weather,selection] = select_typical_72h( ...
    date,T,RH,P,CO2ppm,w)
nHour = numel(T);
startIndex = (1:24:(nHour-72+1))';
nCandidate = numel(startIndex);
feature = zeros(nCandidate,4);
for q = 1:nCandidate
    idx = startIndex(q):(startIndex(q)+71);
    feature(q,:) = [mean(T(idx)),mean(w(idx)), ...
        std(T(idx),1),std(w(idx),1)];
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

score = sum(((feature-center)./spread).^2,2);
[bestScore,bestCandidate] = min(score);
idx = startIndex(bestCandidate):(startIndex(bestCandidate)+71);

weather = struct('Date',date(idx),'T',T(idx),'RH',RH(idx), ...
    'P',P(idx),'CO2ppm',CO2ppm(idx),'w',w(idx));
selection = struct('SelectedStart',date(idx(1)), ...
    'SelectedEnd',date(idx(end)),'CandidateCount',nCandidate, ...
    'Score',bestScore,'MeanT',feature(bestCandidate,1), ...
    'MeanW',feature(bestCandidate,2),'StdT',feature(bestCandidate,3), ...
    'StdW',feature(bestCandidate,4));
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
    humidDuty,dehumidDuty,totalPower_kW] = decode_devices( ...
    Q,mw,F,device)
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

function make_sensitivity_figures(resultDir,gamma,summary, ...
    selectedGamma,trackingTolerance_pct,selectedAtUpperBoundary)
fontName = choose_chinese_font();
if usejava('desktop')
    visibleMode = 'on';
else
    visibleMode = 'off';
end
if selectedAtUpperBoundary
    selectionLabel = '区间内推荐值';
else
    selectionLabel = '推荐增益';
end

fig1 = figure('Visible',visibleMode,'Color','w', ...
    'Position',[100 80 1000 680]);
plot(gamma,summary.NormalizedTrackingRMSE,'-o','LineWidth',1.4);
hold on;
plot(gamma,summary.NormalizedTrackingRMSE_First6h,'-s','LineWidth',1.4);
plot(gamma,summary.NormalizedTrackingRMSE_After6h,'-^','LineWidth',1.4);
xline(selectedGamma,'--',selectionLabel,'LineWidth',1.2);
grid on;
xlabel('协同增益\gamma');
ylabel('综合归一化跟踪RMSE');
title('协同增益对跟踪性能的影响');
legend({'72 h总体','前6 h','6～72 h',selectionLabel}, ...
    'Location','best');
set(gca,'XTick',gamma);
set(findall(fig1,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig1,fullfile(resultDir, ...
    'tracking_metrics_vs_gamma.png'),'Resolution',150);

fig2 = figure('Visible',visibleMode,'Color','w', ...
    'Position',[120 100 1000 680]);
plot(gamma,summary.GraphDisagreementRMSE,'-o','LineWidth',1.4);
hold on;
plot(gamma,summary.SameStagePairRMSE,'-s','LineWidth',1.4);
xline(selectedGamma,'--',selectionLabel,'LineWidth',1.2);
grid on;
xlabel('协同增益\gamma');
ylabel('归一化一致性RMSE');
title('协同增益对仓间协调性的影响');
legend({'图一致性误差','同阶段仓对误差',selectionLabel}, ...
    'Location','best');
set(gca,'XTick',gamma);
set(findall(fig2,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig2,fullfile(resultDir, ...
    'coordination_metrics_vs_gamma.png'),'Resolution',150);

fig3 = figure('Visible',visibleMode,'Color','w', ...
    'Position',[140 120 1050 700]);
feasible = summary.Feasible;
scatter(summary.TrackingChange_pct_vs_gamma0(~feasible), ...
    summary.CoordinationImprovement_pct_vs_gamma0(~feasible), ...
    70,[0.60 0.60 0.60],'filled');
hold on;
scatter(summary.TrackingChange_pct_vs_gamma0(feasible), ...
    summary.CoordinationImprovement_pct_vs_gamma0(feasible), ...
    80,[0.10 0.45 0.75],'filled');
selected = summary.Selected;
scatter(summary.TrackingChange_pct_vs_gamma0(selected), ...
    summary.CoordinationImprovement_pct_vs_gamma0(selected), ...
    150,[0.85 0.20 0.10],'filled','MarkerEdgeColor','k');
xline(trackingTolerance_pct,'--','跟踪容差');
yline(0,':');
for q = 1:numel(gamma)
    text(summary.TrackingChange_pct_vs_gamma0(q)+0.015, ...
        summary.CoordinationImprovement_pct_vs_gamma0(q), ...
        sprintf('\\gamma=%.2f',gamma(q)));
end
grid on;
xlabel('相对\gamma=0的跟踪RMSE变化/%');
ylabel('综合协调改善/%');
title('跟踪性能与仓间协调性的折中关系');
legend({'不可行或独立基准','可行协同增益',selectionLabel}, ...
    'Location','best');
set(findall(fig3,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig3,fullfile(resultDir, ...
    'tracking_coordination_tradeoff.png'),'Resolution',150);

fig4 = figure('Visible',visibleMode,'Color','w', ...
    'Position',[160 140 1050 760]);
subplot(3,1,1);
plot(gamma,summary.SafetyViolationSamples,'-o','LineWidth',1.3);
grid on;
ylabel('越界样本数');
title('安全、饱和与设备能耗随协同增益变化');
subplot(3,1,2);
plot(gamma,summary.SaturationAtControlUpdates,'-o','LineWidth',1.3);
grid on;
ylabel('饱和次数');
subplot(3,1,3);
plot(gamma,summary.DeviceEnergy_kWh,'-o','LineWidth',1.3);
grid on;
xlabel('协同增益\gamma');
ylabel('设备能耗代理/kWh');
set(findall(fig4,'-property','FontName'),'FontName',fontName);
drawnow;
exportgraphics(fig4,fullfile(resultDir, ...
    'constraints_energy_vs_gamma.png'),'Resolution',150);
end

function write_provenance(filePath,weatherFile,selection,gammaGrid, ...
    trackingTolerance_pct,nearMaximumFraction,dt,controlPeriod, ...
    weatherPeriod,tEnd)
Parameter = {'weatherFile';'selectedWeatherStart';'gammaGrid'; ...
    'trackingTolerance';'nearMaximumCoordinationFraction'; ...
    'selectionRule';'dt';'controlPeriod';'weatherPeriod';'tEnd'; ...
    'initialRoom2HumidityRatio';'physicalModel';'controllerH'};
Value = {weatherFile;char(selection.SelectedStart); ...
    sprintf('%.2f ',gammaGrid);sprintf('%.3f',trackingTolerance_pct); ...
    sprintf('%.3f',nearMaximumFraction); ...
    ['tracking/safety/saturation feasible; >=90% maximum coordination ', ...
     'improvement; choose smallest positive gamma']; ...
    sprintf('%g',dt);sprintf('%g',controlPeriod); ...
    sprintf('%g',weatherPeriod);sprintf('%g',tEnd); ...
    '0.0165089';'same as step1 and experiment2'; ...
    'same as experiment1 and experiment2'};
Unit = {'path';'datetime';'-';'%';'-';'-';'s';'s';'s';'s'; ...
    'kg/kgda';'-';'-'};
Category = {'用户提供文件';'派生计算值';'算法设计值'; ...
    '算法设计值';'算法设计值';'算法设计值'; ...
    '用户指定值';'用户指定值';'用户指定值';'用户指定值'; ...
    '用户指定值';'继承模型';'算法设计值'};
SourceNote = {'南宁2023逐时重建气象表'; ...
    '与实验二相同的固定典型窗口规则'; ...
    '仅用于参数敏感性扫描，不是文献推荐值'; ...
    '允许相对独立控制的72 h跟踪增幅'; ...
    '避免为最后少量协调改善选择过大gamma'; ...
    '仿真前固定的可复现选择规则'; ...
    '与实验一和实验二一致';'与实验一和实验二一致'; ...
    '与实验二一致';'与实验二一致'; ...
    '明确排除作废值0.0171547'; ...
    '温度-含湿量-CO2三状态模型'; ...
    '本实验只改变gamma，不改变H_g'};
tableOut = table(Parameter,Value,Unit,Category,SourceNote);
writetable(tableOut,filePath);
end

function fontName = choose_chinese_font
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
