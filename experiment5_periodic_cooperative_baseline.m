function experiment5_periodic_cooperative_baseline
%EXPERIMENT5_PERIODIC_COOPERATIVE_BASELINE 论文仿真实验E5：周期协同控制基准。
%
% 冻结条件：
%   1. 六仓非线性模型，房间阶段为G1、G1、G2、G3、G3、G4；
%   2. gamma=0.35，每仓4台除湿器；
%   3. 物理积分步长60 s，控制与全量通信周期300 s；
%   4. 使用与E4一致的典型72 h动态气象；
%   5. 每个控制时刻六仓均发送最新状态，控制指令零阶保持300 s。
%
% 研究目的：
%   建立后续独立控制、固定事件触发和自适应事件触发实验的统一基准。
%
% 运行前提：本文件应与mushroom_system_parameters.m位于同一目录，且有
%   results/full_precision_matrices/feedback_gains.mat

clc;
close all;
format long g;

%% 1. 文件位置、参数和全精度反馈增益
codeDirectory = fileparts(mfilename('fullpath'));
resultDirectory = fullfile(codeDirectory, ...
    'results','paper_rerun','E5_periodic_cooperative_baseline');
if ~exist(resultDirectory,'dir')
    mkdir(resultDirectory);
end

cfg = mushroom_system_parameters();

assert(cfg.time.dt_s==60, ...
    'E5要求物理积分步长为60 s。');
assert(cfg.time.controlPeriod_s==300, ...
    'E5要求控制周期为300 s。');
assert(cfg.time.weatherPeriod_s==3600, ...
    'E5要求气象更新周期为3600 s。');
assert(cfg.time.mainDuration_s==72*3600, ...
    'E5要求主仿真时长为72 h。');
assert(abs(cfg.control.gamma-0.35)<1e-12, ...
    '唯一参数入口中的协同增益不是冻结值0.35。');
assert(cfg.device.nDehumidifier==4, ...
    'E4确定的每仓除湿器数量应为4台。');

gainFile = fullfile(codeDirectory,'results', ...
    'full_precision_matrices','feedback_gains.mat');
assert(isfile(gainFile), ...
    ['未找到全精度反馈增益文件：\n%s\n', ...
     '请保持results/full_precision_matrices目录结构。'],gainFile);

gainData = load(gainFile,'gain');
assert(isfield(gainData,'gain') && isfield(gainData.gain,'H'), ...
    'feedback_gains.mat中缺少gain.H。');
H = gainData.gain.H;
assert(isequal(size(H),[3 3 4]) && all(isfinite(H(:))), ...
    '反馈增益H的尺寸或数值不正确。');

if ~cfg.device.verified
    warning(['设备规格仍标记为“待核对”。E5可验证控制与容量的仿真可行性，', ...
        '但不能代替厂家数据核对。']);
end

%% 2. 仿真时钟
dt = cfg.time.dt_s;
controlPeriod = cfg.time.controlPeriod_s;
weatherPeriod = cfg.time.weatherPeriod_s;
tEnd = cfg.time.mainDuration_s;

assert(mod(controlPeriod,dt)==0 && mod(weatherPeriod,dt)==0, ...
    '控制周期和气象周期必须为积分步长的整数倍。');

t = (0:dt:tEnd)';
clock = struct;
clock.dt = dt;
clock.controlPeriod = controlPeriod;
clock.weatherPeriod = weatherPeriod;
clock.tEnd = tEnd;
clock.t = t;
clock.nStep = numel(t);
clock.nRoom = cfg.room.nRoom;
clock.nPhysicsPerControl = controlPeriod/dt;

%% 3. 读取并选择与E4一致的典型72 h气象
weatherFile = char(cfg.files.weather);
assert(isfile(weatherFile),'未找到气象文件：\n%s',weatherFile);

raw = readmatrix(weatherFile,'Sheet','第四部分模型输入', ...
    'Range','A2:G8761');
assert(size(raw,1)==8760 && size(raw,2)>=6, ...
    '气象表应含8760行且至少6列，当前为%d x %d。', ...
    size(raw,1),size(raw,2));
assert(all(isfinite(raw(:,1:6)),'all'), ...
    '气象表A:F范围含缺失或非数值数据。');

allDate = datetime(raw(:,1),'ConvertFrom','excel');
allTout = raw(:,2);
allRH = raw(:,3);
allP = raw(:,5);
allCO2ppm = raw(:,6);
allWout = humidity_ratio_from_t_rh_p(allTout,allRH,allP);

[weather,selection] = select_typical_72h( ...
    allDate,allTout,allRH,allP,allCO2ppm,allWout);

weatherHourIndex = min(floor(t/weatherPeriod)+1,72);
outdoor.T = weather.T(weatherHourIndex);
outdoor.RH = weather.RH(weatherHourIndex);
outdoor.w = weather.w(weatherHourIndex);
outdoor.P = weather.P(weatherHourIndex);
outdoor.Cppm = weather.CO2ppm(weatherHourIndex);
outdoor.c = ppm_to_mass(outdoor.Cppm,outdoor.T, ...
    outdoor.P,cfg.physical);

%% 4. 建立模型并运行300 s周期协同控制
model = build_model_from_cfg(cfg,H);
initial.T = cfg.initial.T_C(:)';
initial.w = cfg.initial.w_kg_kgda(:)';
initial.Cppm = cfg.initial.CO2_ppm(:)';

assert(abs(initial.w(2)-0.0165089)<1e-12, ...
    '2号仓初始含湿量不是指定值0.0165089。');

fprintf('\n============================================================\n');
fprintf('E5 六仓周期协同控制基准开始\n');
fprintf('gamma = %.6f，每仓除湿器 = %d台\n', ...
    cfg.control.gamma,cfg.device.nDehumidifier);
fprintf('动态气象 = %s 至 %s\n', ...
    char(selection.SelectedStart),char(selection.SelectedEnd));
fprintf('控制与全量通信周期 = %g s，仿真时长 = %g h\n', ...
    controlPeriod,tEnd/3600);
fprintf('============================================================\n');

baseline = simulate_periodic(initial,model,clock,outdoor);
[baseline,summaryTable,roomSummary] = ...
    analyze_baseline(baseline,model,clock,cfg);

expectedTransmissions = nnz(baseline.isControlUpdate)*clock.nRoom;
AllPeriodicBroadcastsCompleted = ...
    baseline.metric.NetworkTransmissions==expectedTransmissions;
FiniteSimulation = all(isfinite(baseline.T),'all') ...
    && all(isfinite(baseline.w),'all') ...
    && all(isfinite(baseline.Cppm),'all');
ZeroSafetyViolation = baseline.metric.SafetyViolationSamples==0;
NoSaturationAtControlUpdates = ...
    baseline.metric.SaturationAtControlUpdates==0;
NoInterlockConflict = baseline.metric.InterlockConflicts==0;
NonnegativeDehumidificationMargin = ...
    baseline.metric.DehumidificationCapacityMargin_kg_h>=-1e-12;

Pass = FiniteSimulation && AllPeriodicBroadcastsCompleted ...
    && ZeroSafetyViolation && NoSaturationAtControlUpdates ...
    && NoInterlockConflict && NonnegativeDehumidificationMargin;

checkTable = table(FiniteSimulation,AllPeriodicBroadcastsCompleted, ...
    ZeroSafetyViolation,NoSaturationAtControlUpdates, ...
    NoInterlockConflict,NonnegativeDehumidificationMargin,Pass);

%% 5. 保存数值结果
writetable(summaryTable,fullfile(resultDirectory,'baseline_summary.csv'));
writetable(roomSummary,fullfile(resultDirectory,'room_summary.csv'));
writetable(checkTable,fullfile(resultDirectory,'pass_checks.csv'));

write_state_outputs(resultDirectory,baseline,clock,outdoor,cfg,model);
write_device_outputs(resultDirectory,baseline,clock,cfg,model);
write_communication_outputs(resultDirectory,baseline,clock,cfg);

selectedWeatherTable = table((0:71)',weather.Date,weather.T, ...
    weather.RH,weather.P,weather.w,weather.CO2ppm, ...
    'VariableNames',{'HourIndex','DateTime_Beijing','T_out_C', ...
    'RH_out_pct','P_out_Pa','w_out_kg_kgda','CO2_background_ppm'});
writetable(selectedWeatherTable,fullfile(resultDirectory, ...
    'selected_weather_72h.csv'));

write_provenance(fullfile(resultDirectory,'parameter_provenance.csv'), ...
    cfg,weatherFile,gainFile,selection);

save(fullfile(resultDirectory,'E5_periodic_cooperative_baseline.mat'), ...
    'baseline','summaryTable','roomSummary','checkTable','weather', ...
    'selection','model','clock','outdoor','cfg','Pass');

try
    make_figures(resultDirectory,baseline,clock,weather,cfg,model);
catch ME
    warning('数值结果已保存，但图形生成失败：%s',ME.message);
end

%% 6. 命令窗口摘要
fprintf('\n================ E5结果汇总 ================\n');
disp(summaryTable);
disp(roomSummary);

fprintf('检查结果：\n');
fprintf('1. 数值均为有限数：%s\n',logical_text(FiniteSimulation));
fprintf('2. 周期通信完整执行：%s\n', ...
    logical_text(AllPeriodicBroadcastsCompleted));
fprintf('3. 安全越界为0：%s\n',logical_text(ZeroSafetyViolation));
fprintf('4. 控制时刻设备饱和为0：%s\n', ...
    logical_text(NoSaturationAtControlUpdates));
fprintf('5. 设备互锁冲突为0：%s\n', ...
    logical_text(NoInterlockConflict));
fprintf('6. 除湿容量余量非负：%s\n', ...
    logical_text(NonnegativeDehumidificationMargin));

if Pass
    fprintf('\n结论：PASS，E5周期协同控制基准建立完成。\n');
    fprintf('可将本结果作为后续控制方法公平比较的共同基准。\n');
else
    fprintf('\n结论：FAIL，请根据未通过项目检查参数或模型。\n');
end
fprintf('结果目录：%s\n',resultDirectory);
fprintf('=============================================\n');

end


function model = build_model_from_cfg(cfg,H)
model.par = cfg.physical;
model.T0 = cfg.stage.TTarget_C;
model.w0 = cfg.stage.wTarget_kg_kgda;
model.Cppm0 = cfg.stage.CO2Target_ppm;
model.Mda = cfg.stage.Mda_kg;
model.Qbio = cfg.bio.Q_W;
model.mbioW = cfg.bio.mWater_kg_s;
model.mbioC = cfg.bio.mCO2_kg_s;
model.TMin = cfg.safety.TMin_C;
model.TMax = cfg.safety.TMax_C;
model.wMin = cfg.safety.wMin_kg_kgda;
model.wMax = cfg.safety.wMax_kg_kgda;
model.CppmMin = cfg.safety.CO2Min_ppm;
model.CppmMax = cfg.safety.CO2Max_ppm;
model.scaleT = cfg.normalization.scaleT_C;
model.scaleW = cfg.normalization.scaleW_kg_kgda;
model.scaleCppm = cfg.normalization.scaleCO2_ppm;
model.roomStage = cfg.room.stageIndex;
model.A = cfg.network.A;
model.L = cfg.network.L;
model.gamma = cfg.control.gamma;
model.H = H;
model.limits.QMin = cfg.limits.QMin_W;
model.limits.QMax = cfg.limits.QMax_W;
model.limits.mwMin = cfg.limits.mWaterMin_kg_s;
model.limits.mwMax = cfg.limits.mWaterMax_kg_s;
model.limits.FMin = cfg.limits.FMin_m3_s;
model.limits.FMax = cfg.limits.FMax_m3_s;
model.device = cfg.device;
end


function result = simulate_periodic(initial,model,clock,outdoor)
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
        model.par.PIndoorRef_Pa,model.par);
end

Qraw = zeros(nStep,nRoom);
mwRaw = zeros(nStep,nRoom);
Fraw = zeros(nStep,nRoom);
Qcmd = zeros(nStep,nRoom);
mwCmd = zeros(nStep,nRoom);
Fcmd = zeros(nStep,nRoom);
satDetail = false(nStep,nRoom,6);
isControlUpdate = false(nStep,1);
isBroadcast = false(nStep,nRoom);

uRawHold = zeros(nRoom,3);
uHold = zeros(nRoom,3);
satHold = false(nRoom,6);

for k = 1:nStep
    out = struct('T',outdoor.T(k),'w',outdoor.w(k), ...
        'P',outdoor.P(k),'Cppm',outdoor.Cppm(k),'c',outdoor.c(k));

    if mod(k-1,clock.nPhysicsPerControl)==0
        isControlUpdate(k) = true;

        % 周期基准：每个控制时刻六仓全部发送最新状态。
        zCurrent = normalized_state(T(k,:),w(k,:),c(k,:),model);
        isBroadcast(k,:) = true(1,nRoom);

        [uRawHold,uHold,satHold] = ...
            controller_from_current_state(zCurrent,model,out);
    end

    Qraw(k,:) = uRawHold(:,1)';
    mwRaw(k,:) = uRawHold(:,2)';
    Fraw(k,:) = uRawHold(:,3)';
    Qcmd(k,:) = uHold(:,1)';
    mwCmd(k,:) = uHold(:,2)';
    Fcmd(k,:) = uHold(:,3)';
    satDetail(k,:,:) = reshape(satHold,[1,nRoom,6]);

    if k==nStep
        break;
    end

    for i = 1:nRoom
        g = model.roomStage(i);
        xk = [T(k,i);w(k,i);c(k,i)];
        ui = uHold(i,:)';
        bio = [model.Qbio(g);model.mbioW(g);model.mbioC(g)];

        k1 = room_rhs(xk,ui,bio,model.Mda(g),out,model.par);
        k2 = room_rhs(xk+0.5*dt*k1,ui,bio,model.Mda(g),out,model.par);
        k3 = room_rhs(xk+0.5*dt*k2,ui,bio,model.Mda(g),out,model.par);
        k4 = room_rhs(xk+dt*k3,ui,bio,model.Mda(g),out,model.par);
        xNext = xk+dt/6*(k1+2*k2+2*k3+k4);

        T(k+1,i) = xNext(1);
        w(k+1,i) = xNext(2);
        c(k+1,i) = max(xNext(3),0);
    end
end

result = struct('T',T,'w',w,'c',c,'Qraw',Qraw, ...
    'mwRaw',mwRaw,'Fraw',Fraw,'Qcmd',Qcmd,'mwCmd',mwCmd, ...
    'Fcmd',Fcmd,'satDetail',satDetail, ...
    'isControlUpdate',isControlUpdate,'isBroadcast',isBroadcast);
end


function z = normalized_state(Trow,wrow,crow,model)
nRoom = numel(model.roomStage);
z = zeros(nRoom,3);
for i = 1:nRoom
    g = model.roomStage(i);
    CppmNow = mass_to_ppm(crow(i),Trow(i), ...
        model.par.PIndoorRef_Pa,model.par);
    z(i,:) = [(Trow(i)-model.T0(g))/model.scaleT(g), ...
              (wrow(i)-model.w0(g))/model.scaleW(g), ...
              (CppmNow-model.Cppm0(g))/model.scaleCppm(g)];
end
end


function [uRaw,uApplied,sat] = ...
    controller_from_current_state(zCurrent,model,out)
% eta包含本仓跟踪偏差与环形网络邻仓偏差的一致性校正。
eta = zCurrent+model.gamma*(model.L*zCurrent);
uFF = equilibrium_feedforward(model,out);
nRoom = numel(model.roomStage);
uRaw = zeros(nRoom,3);
for i = 1:nRoom
    g = model.roomStage(i);
    uRaw(i,:) = uFF(g,:)+(model.H(:,:,g)*eta(i,:)')';
end

uApplied = uRaw;
uApplied(:,1) = min(max(uRaw(:,1),model.limits.QMin), ...
    model.limits.QMax);
uApplied(:,2) = min(max(uRaw(:,2),model.limits.mwMin), ...
    model.limits.mwMax);
uApplied(:,3) = min(max(uRaw(:,3),model.limits.FMin), ...
    model.limits.FMax);

sat = [uRaw(:,1)<model.limits.QMin, ...
       uRaw(:,1)>model.limits.QMax, ...
       uRaw(:,2)<model.limits.mwMin, ...
       uRaw(:,2)>model.limits.mwMax, ...
       uRaw(:,3)<model.limits.FMin, ...
       uRaw(:,3)>model.limits.FMax];
end


function uFF = equilibrium_feedforward(model,out)
nStage = numel(model.T0);
cTarget = zeros(nStage,1);
for g = 1:nStage
    cTarget(g) = ppm_to_mass(model.Cppm0(g),model.T0(g), ...
        model.par.PIndoorRef_Pa,model.par);
end
F0 = model.mbioC./(cTarget-out.c);
mw0 = model.par.rhoAir_kg_m3.*F0.*(model.w0-out.w) ...
    -model.mbioW;
Q0 = (model.par.UA_W_K+model.par.rhoAir_kg_m3.* ...
    model.par.cpAir_J_kgK.*F0).*(model.T0-out.T)-model.Qbio;
uFF = [Q0,mw0,F0];
end


function dx = room_rhs(x,u,bio,Mda,out,par)
T = x(1);
w = x(2);
c = x(3);
Q = u(1);
mw = u(2);
F = u(3);
dT = (par.UA_W_K*(out.T-T) ...
    +par.rhoAir_kg_m3*par.cpAir_J_kgK*F*(out.T-T) ...
    +Q+bio(1))/par.CT_J_K;
dw = (par.rhoAir_kg_m3*F*(out.w-w)+mw+bio(2))/Mda;
dc = (F*(out.c-c)+bio(3))/par.V_m3;
dx = [dT;dw;dc];
end


function [result,summaryTable,roomSummary] = ...
    analyze_baseline(result,model,clock,cfg)
nStep = clock.nStep;
nRoom = clock.nRoom;

Cppm = mass_to_ppm(result.c,result.T, ...
    model.par.PIndoorRef_Pa,model.par);

targetT = model.T0(model.roomStage)';
targetW = model.w0(model.roomStage)';
targetC = model.Cppm0(model.roomStage)';
safeTMin = model.TMin(model.roomStage)';
safeTMax = model.TMax(model.roomStage)';
safeWMin = model.wMin(model.roomStage)';
safeWMax = model.wMax(model.roomStage)';
safeCMin = model.CppmMin(model.roomStage)';
safeCMax = model.CppmMax(model.roomStage)';
normT = model.scaleT(model.roomStage)';
normW = model.scaleW(model.roomStage)';
normC = model.scaleCppm(model.roomStage)';

dT = result.T-targetT;
dW = result.w-targetW;
dC = Cppm-targetC;
zT = dT./normT;
zW = dW./normW;
zC = dC./normC;
z = cat(3,zT,zW,zC);
errorNorm = sqrt(sum(z.^2,3));

violT = result.T<safeTMin | result.T>safeTMax;
violW = result.w<safeWMin | result.w>safeWMax;
violC = Cppm<safeCMin | Cppm>safeCMax;
anyViolation = violT|violW|violC;

[cool,heat,humid,dehumid,~,~,~,~,~,~,~,devicePower] = ...
    decode_devices(result.Qcmd,result.mwCmd,result.Fcmd,model.device);
interlockConflicts = nnz(cool>1e-12 & heat>1e-12) ...
    +nnz(humid>1e-12 & dehumid>1e-12);

controlIndex = find(result.isControlUpdate);
satAtUpdates = result.satDetail(controlIndex,:,:);
satQAtUpdates = satAtUpdates(:,:,1)|satAtUpdates(:,:,2);
satWAtUpdates = satAtUpdates(:,:,3)|satAtUpdates(:,:,4);
satFAtUpdates = satAtUpdates(:,:,5)|satAtUpdates(:,:,6);

nBroadcast = nnz(result.isBroadcast);
nReceive = sum(model.A,'all')/nRoom*nBroadcast;
communicationEnergy_J = (0.366*nBroadcast+0.310*nReceive)/1000;
deviceEnergy_kWh = sum(devicePower(1:end-1,:),'all') ...
    *clock.dt/3600;

graphResidual = zeros(size(z));
for k = 1:nStep
    graphResidual(k,:,:) = model.L*squeeze(z(k,:,:));
end
sameStageDifference = cat(2,z(:,1,:)-z(:,2,:), ...
    z(:,4,:)-z(:,5,:));

capacity_kg_h = -model.limits.mwMin*3600;
peakRequestedDeh_kg_h = max(max(-result.mwRaw*3600,0),[],'all');
peakAppliedDeh_kg_h = max(max(-result.mwCmd*3600,0),[],'all');
capacityMargin_kg_h = capacity_kg_h-peakRequestedDeh_kg_h;

metric = struct;
metric.NormalizedTrackingRMSE = sqrt(mean(z.^2,'all'));
metric.TemperatureRMSE_C = sqrt(mean(dT.^2,'all'));
metric.HumidityRatioRMSE = sqrt(mean(dW.^2,'all'));
metric.CO2RMSE_ppm = sqrt(mean(dC.^2,'all'));
metric.MaximumTemperatureError_C = max(abs(dT),[],'all');
metric.MaximumHumidityRatioError = max(abs(dW),[],'all');
metric.MaximumCO2Error_ppm = max(abs(dC),[],'all');
metric.GraphDisagreementRMSE = sqrt(mean(graphResidual.^2,'all'));
metric.SameStagePairRMSE = sqrt(mean(sameStageDifference.^2,'all'));
metric.SafetyViolationSamples = nnz(anyViolation);
metric.TemperatureViolationSamples = nnz(violT);
metric.HumidityViolationSamples = nnz(violW);
metric.CO2ViolationSamples = nnz(violC);
metric.ThermalSaturationAtControlUpdates = nnz(satQAtUpdates);
metric.WaterSaturationAtControlUpdates = nnz(satWAtUpdates);
metric.VentilationSaturationAtControlUpdates = nnz(satFAtUpdates);
metric.SaturationAtControlUpdates = nnz(satQAtUpdates) ...
    +nnz(satWAtUpdates)+nnz(satFAtUpdates);
metric.InterlockConflicts = interlockConflicts;
metric.ControlUpdates = numel(controlIndex);
metric.NetworkTransmissions = nBroadcast;
metric.NetworkReceptions = nReceive;
metric.ActiveCommunicationEnergy_J = communicationEnergy_J;
metric.DeviceEnergy_kWh = deviceEnergy_kWh;
metric.DehumidificationCapacity_kg_h = capacity_kg_h;
metric.PeakRequestedDehumidification_kg_h = peakRequestedDeh_kg_h;
metric.PeakAppliedDehumidification_kg_h = peakAppliedDeh_kg_h;
metric.DehumidificationCapacityMargin_kg_h = capacityMargin_kg_h;

Condition = {'300 s周期协同控制'};
Gamma = model.gamma;
DehumidifierCountPerRoom = model.device.nDehumidifier;
NormalizedTrackingRMSE = metric.NormalizedTrackingRMSE;
TemperatureRMSE_C = metric.TemperatureRMSE_C;
HumidityRatioRMSE = metric.HumidityRatioRMSE;
CO2RMSE_ppm = metric.CO2RMSE_ppm;
GraphDisagreementRMSE = metric.GraphDisagreementRMSE;
SameStagePairRMSE = metric.SameStagePairRMSE;
SafetyViolationSamples = metric.SafetyViolationSamples;
SaturationAtControlUpdates = metric.SaturationAtControlUpdates;
InterlockConflicts = metric.InterlockConflicts;
NetworkTransmissions = metric.NetworkTransmissions;
NetworkReceptions = metric.NetworkReceptions;
ActiveCommunicationEnergy_J = metric.ActiveCommunicationEnergy_J;
DeviceEnergy_kWh = metric.DeviceEnergy_kWh;
DehumidificationCapacity_kg_h = metric.DehumidificationCapacity_kg_h;
PeakRequestedDehumidification_kg_h = ...
    metric.PeakRequestedDehumidification_kg_h;
DehumidificationCapacityMargin_kg_h = ...
    metric.DehumidificationCapacityMargin_kg_h;

summaryTable = table(Condition,Gamma,DehumidifierCountPerRoom, ...
    NormalizedTrackingRMSE,TemperatureRMSE_C,HumidityRatioRMSE, ...
    CO2RMSE_ppm,GraphDisagreementRMSE,SameStagePairRMSE, ...
    SafetyViolationSamples,SaturationAtControlUpdates, ...
    InterlockConflicts,NetworkTransmissions,NetworkReceptions, ...
    ActiveCommunicationEnergy_J,DeviceEnergy_kWh, ...
    DehumidificationCapacity_kg_h, ...
    PeakRequestedDehumidification_kg_h, ...
    DehumidificationCapacityMargin_kg_h);

Room = cfg.room.names;
Stage = cfg.stage.names(cfg.room.stageIndex);
TemperatureRMSE_C = sqrt(mean(dT.^2,1))';
HumidityRatioRMSE = sqrt(mean(dW.^2,1))';
CO2RMSE_ppm = sqrt(mean(dC.^2,1))';
MaximumNormalizedError = max(errorNorm,[],1)';
SafetyViolationSamples = sum(anyViolation,1)';
SaturationAtControlUpdates = squeeze(sum(satQAtUpdates| ...
    satWAtUpdates|satFAtUpdates,1))';
RecoveryTime_h = recovery_time(z,clock.dt,0.05,3600);

roomSummary = table(Room,Stage,TemperatureRMSE_C, ...
    HumidityRatioRMSE,CO2RMSE_ppm,MaximumNormalizedError, ...
    RecoveryTime_h,SafetyViolationSamples,SaturationAtControlUpdates);

result.Cppm = Cppm;
result.dT = dT;
result.dW = dW;
result.dC = dC;
result.z = z;
result.errorNorm = errorNorm;
result.graphResidual = graphResidual;
result.violT = violT;
result.violW = violW;
result.violC = violC;
result.anyViolation = anyViolation;
result.devicePower_kW = devicePower;
result.metric = metric;
end


function recoveryTime_h = recovery_time(z,dt,tolerance,dwell_s)
nRoom = size(z,2);
dwellSample = round(dwell_s/dt)+1;
recoveryTime_h = NaN(nRoom,1);
for i = 1:nRoom
    within = squeeze(all(abs(z(:,i,:))<=tolerance,3));
    if numel(within)>=dwellSample
        count = conv(double(within),ones(dwellSample,1),'valid');
        first = find(count==dwellSample,1,'first');
        if ~isempty(first)
            recoveryTime_h(i) = (first-1)*dt/3600;
        end
    end
end
end


function write_state_outputs(resultDirectory,result,clock,outdoor,cfg,model)
nStep = clock.nStep;
nRoom = clock.nRoom;
time_s = repmat(clock.t,nRoom,1);
time_h = time_s/3600;
roomIndex = kron((1:nRoom)',ones(nStep,1));
Room = cfg.room.names(roomIndex);
Stage = cfg.stage.names(cfg.room.stageIndex(roomIndex));

TargetT_C = model.T0(model.roomStage(roomIndex));
TargetW_kg_kgda = model.w0(model.roomStage(roomIndex));
TargetCO2_ppm = model.Cppm0(model.roomStage(roomIndex));

stateTable = table(time_s,time_h,Room,Stage,result.T(:), ...
    result.w(:),result.Cppm(:),TargetT_C,TargetW_kg_kgda, ...
    TargetCO2_ppm,result.dT(:),result.dW(:),result.dC(:), ...
    reshape(result.z(:,:,1),[],1),reshape(result.z(:,:,2),[],1), ...
    reshape(result.z(:,:,3),[],1),result.errorNorm(:), ...
    result.violT(:),result.violW(:),result.violC(:), ...
    result.anyViolation(:),repmat(outdoor.T,nRoom,1), ...
    repmat(outdoor.w,nRoom,1),repmat(outdoor.Cppm,nRoom,1), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','T_C', ...
    'w_kg_kgda','CO2_ppm','TargetT_C','TargetW_kg_kgda', ...
    'TargetCO2_ppm','DeltaT_C','DeltaW_kg_kgda','DeltaCO2_ppm', ...
    'NormalizedTError','NormalizedWError','NormalizedCO2Error', ...
    'NormalizedErrorNorm','TemperatureViolation', ...
    'HumidityViolation','CO2Violation','AnySafetyViolation', ...
    'Tout_C','Wout_kg_kgda','CO2out_ppm'});
writetable(stateTable,fullfile(resultDirectory,'state_trajectories.csv'));

writetable(stateTable(stateTable.AnySafetyViolation,:), ...
    fullfile(resultDirectory,'constraint_violations.csv'));
end


function write_device_outputs(resultDirectory,result,clock,cfg,model)
nStep = clock.nStep;
nRoom = clock.nRoom;
time_s = repmat(clock.t,nRoom,1);
time_h = time_s/3600;
roomIndex = kron((1:nRoom)',ones(nStep,1));
Room = cfg.room.names(roomIndex);
Stage = cfg.stage.names(cfg.room.stageIndex(roomIndex));

[cool,heat,humid,dehumid,fanLow,fanHigh,fanOff,coolDuty, ...
    heatDuty,humidDuty,dehumidDuty,power] = ...
    decode_devices(result.Qcmd,result.mwCmd,result.Fcmd,model.device);

deviceTable = table(time_s,time_h,Room,Stage,result.Qraw(:), ...
    result.Qcmd(:),cool(:),heat(:),coolDuty(:),heatDuty(:), ...
    result.mwRaw(:),result.mwCmd(:),humid(:),dehumid(:), ...
    humidDuty(:),dehumidDuty(:),result.Fraw(:),result.Fcmd(:), ...
    fanOff(:),fanLow(:),fanHigh(:),power(:), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','QRaw_W', ...
    'QApplied_W','Cooling_kW','Heating_kW','CoolingDuty', ...
    'HeatingDuty','WaterRaw_kg_s','WaterApplied_kg_s', ...
    'Humidification_kg_h','Dehumidification_kg_h', ...
    'HumidificationDuty','DehumidificationDuty','FlowRaw_m3_s', ...
    'FlowApplied_m3_s','FanOffDuty','FanLowDuty','FanHighDuty', ...
    'EstimatedPower_kW'});
writetable(deviceTable,fullfile(resultDirectory,'device_commands.csv'));
end


function write_communication_outputs(resultDirectory,result,clock,cfg)
updateIndex = find(result.isControlUpdate);
nUpdate = numel(updateIndex);
nRoom = clock.nRoom;
time_s = repmat(clock.t(updateIndex),nRoom,1);
time_h = time_s/3600;
roomIndex = kron((1:nRoom)',ones(nUpdate,1));
Room = cfg.room.names(roomIndex);
Broadcast = reshape(result.isBroadcast(updateIndex,:),[],1);

communicationTable = table(time_s,time_h,Room,Broadcast, ...
    'VariableNames',{'Time_s','Time_h','Room','Broadcast'});
writetable(communicationTable,fullfile(resultDirectory, ...
    'communication_at_control_updates.csv'));
end


function write_provenance(filePath,cfg,weatherFile,gainFile,selection)
Parameter = {'gamma';'dehumidifierCountPerRoom'; ...
    'dehumidifierCapacityPerUnit';'physicsStep';'controlPeriod'; ...
    'weatherPeriod';'duration';'weatherFile';'weatherStart'; ...
    'weatherEnd';'feedbackGainFile';'deviceVerificationStatus'};
Value = {sprintf('%.15g',cfg.control.gamma); ...
    sprintf('%d',cfg.device.nDehumidifier); ...
    sprintf('%.15g',cfg.device.dehPerUnit_kg_h); ...
    sprintf('%g',cfg.time.dt_s);sprintf('%g',cfg.time.controlPeriod_s); ...
    sprintf('%g',cfg.time.weatherPeriod_s); ...
    sprintf('%g',cfg.time.mainDuration_s);weatherFile; ...
    char(selection.SelectedStart);char(selection.SelectedEnd); ...
    gainFile;cfg.device.status};
Unit = {'-';'台/仓';'kg/h';'s';'s';'s';'s';'path'; ...
    '北京时间';'北京时间';'path';'-'};
Category = {'算法设计值';'E4容量验证值';'继承设备工况值'; ...
    '算法设计值';'算法设计值';'算法设计值';'算法设计值'; ...
    '用户提供数据';'派生选择值';'派生选择值'; ...
    'E2/E3全精度结果';'待核对标记'};
writetable(table(Parameter,Value,Unit,Category),filePath);
end


function make_figures(resultDirectory,result,clock,weather,cfg,model)
fontName = choose_chinese_font();
colors = lines(clock.nRoom);
timeHour = clock.t/3600;

fig1 = figure('Color','w','Position',[80 60 1180 860]);
subplot(3,1,1);
plot_state_with_targets(timeHour,result.T, ...
    model.T0(model.roomStage)',colors);
ylabel('温度/℃'); title('六仓温度轨迹与阶段目标');
legend(cfg.room.names,'Location','eastoutside');

subplot(3,1,2);
plot_state_with_targets(timeHour,result.w, ...
    model.w0(model.roomStage)',colors);
ylabel('含湿量/(kg/kg干空气)'); title('六仓含湿量轨迹与阶段目标');

subplot(3,1,3);
plot_state_with_targets(timeHour,result.Cppm, ...
    model.Cppm0(model.roomStage)',colors);
ylabel('CO2/ppm'); xlabel('时间/h');
title('六仓CO2轨迹与阶段目标');
set(findall(fig1,'-property','FontName'),'FontName',fontName);
exportgraphics(fig1,fullfile(resultDirectory, ...
    '图1_周期协同控制六仓状态轨迹.png'),'Resolution',300);

fig2 = figure('Color','w','Position',[100 80 1080 650]);
plot(timeHour,result.errorNorm,'LineWidth',1.1); hold on;
yline(0.05,'--k','5%目标附近判据','LineWidth',1.1);
grid on; box on;
xlabel('时间/h'); ylabel('综合归一化误差范数');
title('300 s周期协同控制下的六仓归一化跟踪误差');
legend(cfg.room.names,'Location','eastoutside');
set(findall(fig2,'-property','FontName'),'FontName',fontName);
exportgraphics(fig2,fullfile(resultDirectory, ...
    '图2_周期协同控制归一化跟踪误差.png'),'Resolution',300);

fig3 = figure('Color','w','Position',[120 60 1160 850]);
subplot(3,1,1);
plot(timeHour,result.Qcmd/1000,'LineWidth',0.9); hold on;
yline(model.limits.QMin/1000,'--r','制冷能力下限');
yline(model.limits.QMax/1000,'--r','制热能力上限');
grid on; box on; ylabel('热量指令/kW');
title('热量调节指令（负值制冷，正值加热）');

subplot(3,1,2);
plot(timeHour,result.mwCmd*3600,'LineWidth',0.9); hold on;
yline(model.limits.mwMin*3600,'--r','除湿能力下限');
yline(model.limits.mwMax*3600,'--r','加湿能力上限');
grid on; box on; ylabel('水分指令/(kg/h)');
title('水分调节指令（负值除湿，正值加湿）');

subplot(3,1,3);
plot(timeHour,result.Fcmd*3600,'LineWidth',0.9); hold on;
yline(model.limits.FMax*3600,'--r','通风能力上限');
grid on; box on; ylabel('通风量/(m^3/h)'); xlabel('时间/h');
title('新风通风指令');
legend(cfg.room.names,'Location','eastoutside');
set(findall(fig3,'-property','FontName'),'FontName',fontName);
exportgraphics(fig3,fullfile(resultDirectory, ...
    '图3_周期协同控制设备指令.png'),'Resolution',300);

fig4 = figure('Color','w','Position',[140 50 1120 900]);
weatherTime = (0:71)';
subplot(4,1,1); stairs(weatherTime,weather.T,'LineWidth',1.2);
grid on; ylabel('℃'); title('典型72 h室外温度');
subplot(4,1,2); stairs(weatherTime,weather.RH,'LineWidth',1.2);
grid on; ylabel('%'); title('室外相对湿度');
subplot(4,1,3); stairs(weatherTime,weather.w,'LineWidth',1.2);
grid on; ylabel('kg/kg干空气'); title('室外含湿量');
subplot(4,1,4); stairs(weatherTime,weather.CO2ppm,'LineWidth',1.2);
grid on; ylabel('ppm'); xlabel('时间/h'); title('室外CO2背景');
set(findall(fig4,'-property','FontName'),'FontName',fontName);
exportgraphics(fig4,fullfile(resultDirectory, ...
    '图4_周期基准所用典型72h气象.png'),'Resolution',300);
end


function plot_state_with_targets(timeHour,state,target,colors)
hold on;
for i = 1:size(state,2)
    plot(timeHour,state(:,i),'Color',colors(i,:), ...
        'LineWidth',1.05);
    yline(target(i),'--','Color',colors(i,:), ...
        'LineWidth',0.7,'HandleVisibility','off');
end
grid on;
box on;
end


function [weather,selection] = select_typical_72h(date,T,RH,P,CO2ppm,w)
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
        center(j) = local_quantile(feature(:,j),0.5);
    end
    spread(j) = local_quantile(feature(:,j),0.75) ...
        -local_quantile(feature(:,j),0.25);
end
assert(all(spread>0),'典型气象窗口IQR标准化失败。');
score = sum(((feature-center)./spread).^2,2);
[minimumScore,bestCandidate] = min(score);
idx = startIndex(bestCandidate):(startIndex(bestCandidate)+71);

weather = struct('Date',date(idx),'T',T(idx),'RH',RH(idx), ...
    'P',P(idx),'CO2ppm',CO2ppm(idx),'w',w(idx));
selection = struct('SelectedStart',date(idx(1)), ...
    'SelectedEnd',date(idx(end)), ...
    'CandidateCount',nCandidate,'Score',minimumScore);
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
    humidDuty,dehumidDuty,totalPower_kW] = ...
    decode_devices(Q,mw,F,device)
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
fanLowDuty(low) = flow(low)/device.fanLow_m3_h;
high = ~low;
fanHighDuty(high) = (flow(high)-device.fanLow_m3_h) ...
    /(device.fanHigh_m3_h-device.fanLow_m3_h);
fanLowDuty(high) = 1-fanHighDuty(high);
fanOffDuty = max(0,1-fanLowDuty-fanHighDuty);

totalPower_kW = device.coolPowerRated_kW*coolDuty ...
    +device.heatPowerRated_kW*heatDuty ...
    +device.humPowerRated_kW*humidDuty ...
    +(device.nDehumidifier*device.dehPowerPerUnit_kW)*dehumidDuty ...
    +device.fanLowPower_kW*fanLowDuty ...
    +device.fanHighPower_kW*fanHighDuty;
end


function c = ppm_to_mass(Cppm,TdegC,pressurePa,physical)
c = Cppm.*1e-6.*pressurePa.*physical.MCO2_kg_mol ...
    ./(physical.R_J_molK.*(TdegC+273.15));
end


function Cppm = mass_to_ppm(c,TdegC,pressurePa,physical)
Cppm = c.*physical.R_J_molK.*(TdegC+273.15) ...
    ./(pressurePa.*physical.MCO2_kg_mol).*1e6;
end


function fontName = choose_chinese_font
fontName = get(groot,'defaultAxesFontName');
availableFonts = listfonts;
preferred = {'PingFang SC','Microsoft YaHei','SimHei', ...
    'Noto Sans CJK SC','Heiti SC','Songti SC','Arial Unicode MS'};
for k = 1:numel(preferred)
    if any(strcmpi(availableFonts,preferred{k}))
        fontName = preferred{k};
        return;
    end
end
end


function textOut = logical_text(value)
if value
    textOut = '通过';
else
    textOut = '未通过';
end
end
