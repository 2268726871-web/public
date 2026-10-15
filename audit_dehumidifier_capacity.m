function audit_dehumidifier_capacity
%AUDIT_DEHUMIDIFIER_CAPACITY 论文仿真实验E4：设备容量可行性核查。
%
% 功能：
%   1. 继承唯一参数入口、全精度反馈增益和非线性六仓模型；
%   2. 固定动态气象、初始状态、gamma=0.35及全部控制参数；
%   3. 只改变每仓除湿器台数，依次扫描3、4、5台；
%   4. 比较安全越界、除湿饱和、跟踪误差和设备能耗；
%   5. 自动选择“零越界且无除湿饱和”的最小设备配置；
%   6. 输出容量扫描表、逐仓越界表、状态曲线和容量需求图。
%
% 运行前提：
%   本文件应放在“蘑菇仓/代码库”目录内，与下列文件同目录：
%       mushroom_system_parameters.m
%       results/full_precision_matrices/feedback_gains.mat
%
% 注意：本程序不会自动修改mushroom_system_parameters.m。
%
% 数值性质：
%   [I] 继承值；[D] 派生计算值；[A] 算法设计值；[P] 待核对值。
%
% 候选台数属于算法设计值；单台除湿能力及功率仍需设备资料核对。

clc;
close all;
format long g;

%% 1. 文件位置和唯一参数入口
codeDirectory = fileparts(mfilename('fullpath'));
resultDirectory = fullfile(codeDirectory, ...
    'results','paper_rerun','E4_device_capacity');
if ~exist(resultDirectory,'dir')
    mkdir(resultDirectory);
end

cfg = mushroom_system_parameters();

matrixDirectory = fullfile(codeDirectory, ...
    'results','full_precision_matrices');
gainFile = fullfile(matrixDirectory,'feedback_gains.mat');

assert(isfile(gainFile), ...
    ['未找到全精度反馈增益文件：\n%s\n', ...
     '请先运行build_feedback_gains。'],gainFile);

gainData = load(gainFile,'gain');

assert(isfield(gainData,'gain') && isfield(gainData.gain,'H'), ...
    '反馈增益文件中缺少gain.H。');

H = gainData.gain.H;                         % [D]

assert(isequal(size(H),[3 3 4]) && all(isfinite(H(:))), ...
    '全精度反馈增益H的尺寸或数值不正确。');
assert(abs(cfg.control.gamma-0.35)<1e-12, ...
    '唯一参数入口中的协同增益不是本轮冻结值0.35。');

%% 2. 除湿设备容量候选
% [A] 只比较3、4、5台；不改变单台额定除湿能力。
unitCandidates = (3:5)';

% [A] 周期通信不使用事件阈值，该容差只用于数值接口兼容。
eventNumericalTolerance = 1e-12;

assert(cfg.device.dehPerUnit_kg_h>0, ...
    '单台除湿能力必须为正数。');

%% 3. 时间设置
dt = cfg.time.dt_s;
controlPeriod = cfg.time.controlPeriod_s;
eventCheckPeriod = cfg.time.eventCheckPeriod_s;
weatherPeriod = cfg.time.weatherPeriod_s;
tEnd = cfg.time.mainDuration_s;

assert(controlPeriod==eventCheckPeriod, ...
    '本实验要求控制周期和事件检查周期均为300 s。');
assert(mod(controlPeriod,dt)==0 && mod(weatherPeriod,dt)==0, ...
    '控制周期和气象周期必须为积分步长的整数倍。');

t = (0:dt:tEnd)';
clock = struct;
clock.dt = dt;
clock.controlPeriod = controlPeriod;
clock.eventCheckPeriod = eventCheckPeriod;
clock.weatherPeriod = weatherPeriod;
clock.tEnd = tEnd;
clock.t = t;
clock.nStep = numel(t);
clock.nRoom = cfg.room.nRoom;
clock.nPhysicsPerCheck = eventCheckPeriod/dt;

%% 4. 读取并固定典型72 h动态气象
weatherFile = char(cfg.files.weather);
assert(isfile(weatherFile), ...
    '未找到气象文件：\n%s',weatherFile);

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
allWout = humidity_ratio_from_t_rh_p( ...
    allTout,allRH,allP);

[weather,selection] = select_typical_72h( ...
    allDate,allTout,allRH,allP,allCO2ppm,allWout);

weatherHourIndex = min(floor(t/weatherPeriod)+1,72);
outdoor.T = weather.T(weatherHourIndex);
outdoor.w = weather.w(weatherHourIndex);
outdoor.P = weather.P(weatherHourIndex);
outdoor.Cppm = weather.CO2ppm(weatherHourIndex);
outdoor.c = ppm_to_mass(outdoor.Cppm,outdoor.T, ...
    outdoor.P,cfg.physical);

%% 5. 由唯一参数入口建立非线性仿真模型
model = build_model_from_cfg(cfg,H);

initial.T = cfg.initial.T_C(:)';
initial.w = cfg.initial.w_kg_kgda(:)';
initial.Cppm = cfg.initial.CO2_ppm(:)';

assert(abs(initial.w(2)-0.0165089)<1e-12, ...
    '2号仓初始含湿量不是指定的0.0165089。');

if ~cfg.device.verified
    warning(['当前设备规格仍标记为“待核对”。', ...
        '本程序给出的是容量可行性仿真结果，不能替代厂商资料。']);
end

%% 6. 依次运行3、4、5台除湿器的周期通信仿真
nCandidate = numel(unitCandidates);
capacityCases = cell(nCandidate,1);

UnitCount = unitCandidates;
PerUnitCapacity_kg_h = repmat( ...
    cfg.device.dehPerUnit_kg_h,nCandidate,1);
TotalCapacity_kg_h = UnitCount.*PerUnitCapacity_kg_h;

PeakRequestedDehumidification_kg_h = zeros(nCandidate,1);
CapacityMargin_kg_h = zeros(nCandidate,1);
NormalizedTrackingRMSE = zeros(nCandidate,1);
TemperatureRMSE_C = zeros(nCandidate,1);
HumidityRatioRMSE = zeros(nCandidate,1);
CO2RMSE_ppm = zeros(nCandidate,1);
SafetyViolationSamples = zeros(nCandidate,1);
TemperatureViolationSamples = zeros(nCandidate,1);
HumidityViolationSamples = zeros(nCandidate,1);
CO2ViolationSamples = zeros(nCandidate,1);
DehumidificationSaturationAtChecks = zeros(nCandidate,1);
AllDeviceSaturationAtChecks = zeros(nCandidate,1);
InterlockConflicts = zeros(nCandidate,1);
MaximumHumidityRatio = zeros(nCandidate,1);
MaximumHumidityExcess = zeros(nCandidate,1);
DeviceEnergy_kWh = zeros(nCandidate,1);

roomViolationCount = zeros(nCandidate,clock.nRoom);
roomHumidityViolationCount = zeros(nCandidate,clock.nRoom);

fprintf('\n====================================================\n');
fprintf('除湿设备容量可行性核查开始\n');
fprintf('动态气象 = %s 至 %s\n', ...
    char(selection.SelectedStart),char(selection.SelectedEnd));
fprintf('单台除湿能力 = %.6f kg/h\n', ...
    cfg.device.dehPerUnit_kg_h);
fprintf('候选台数 = %s（每仓）\n',mat2str(unitCandidates'));
fprintf('====================================================\n');

for q = 1:nCandidate
    nUnit = unitCandidates(q);
    modelCase = model;

    % 只修改每仓设备台数以及由此派生的除湿输入下限。
    modelCase.device.nDehumidifier = nUnit;
    modelCase.limits.mwMin = ...
        -(nUnit*modelCase.device.dehPerUnit_kg_h)/3600;

    fprintf('\n正在运行每仓%d台除湿器：总能力%.4f kg/h……\n', ...
        nUnit,-modelCase.limits.mwMin*3600);

    oneCase = simulate_case('Periodic',NaN,initial,modelCase, ...
        clock,outdoor,eventNumericalTolerance);
    oneCase = analyze_case(oneCase,modelCase,clock);
    capacityCases{q} = oneCase;

    m = oneCase.metric;
    checkIndex = find(oneCase.isEventCheck);
    requestedDehumidification = max(-oneCase.mwRaw*3600,0);

    PeakRequestedDehumidification_kg_h(q) = ...
        max(requestedDehumidification,[],'all');
    CapacityMargin_kg_h(q) = TotalCapacity_kg_h(q) ...
        - PeakRequestedDehumidification_kg_h(q);

    NormalizedTrackingRMSE(q) = m.NormalizedTrackingRMSE;
    TemperatureRMSE_C(q) = m.TemperatureRMSE_C;
    HumidityRatioRMSE(q) = m.HumidityRatioRMSE;
    CO2RMSE_ppm(q) = m.CO2RMSE_ppm;
    SafetyViolationSamples(q) = m.SafetyViolationSamples;
    TemperatureViolationSamples(q) = m.TemperatureViolationSamples;
    HumidityViolationSamples(q) = m.HumidityViolationSamples;
    CO2ViolationSamples(q) = m.CO2ViolationSamples;
    AllDeviceSaturationAtChecks(q) = m.SaturationAtChecks;
    InterlockConflicts(q) = m.InterlockConflicts;
    DeviceEnergy_kWh(q) = m.DeviceEnergy_kWh;

    dehSaturation = oneCase.mwRaw(checkIndex,:) ...
        < modelCase.limits.mwMin-1e-12;
    DehumidificationSaturationAtChecks(q) = nnz(dehSaturation);

    safeWMaxRoom = modelCase.wMax(modelCase.roomStage)';
    humidityExcess = max(oneCase.w-safeWMaxRoom,0);
    MaximumHumidityRatio(q) = max(oneCase.w,[],'all');
    MaximumHumidityExcess(q) = max(humidityExcess,[],'all');

    for i = 1:clock.nRoom
        roomViolationCount(q,i) = nnz(oneCase.anyViolation(:,i));
        roomHumidityViolationCount(q,i) = nnz(oneCase.violW(:,i));
    end
end

%% 7. 建立容量可行性表并选择最小合格配置
ZeroSafetyViolation = SafetyViolationSamples==0;
NoDehumidificationSaturation = ...
    DehumidificationSaturationAtChecks==0;
NoInterlockConflict = InterlockConflicts==0;

SafetyFeasible = ZeroSafetyViolation & NoInterlockConflict;
CapacityFeasible = SafetyFeasible ...
    & NoDehumidificationSaturation ...
    & CapacityMargin_kg_h>=-1e-12;

capacityTable = table(UnitCount,PerUnitCapacity_kg_h, ...
    TotalCapacity_kg_h,PeakRequestedDehumidification_kg_h, ...
    CapacityMargin_kg_h,NormalizedTrackingRMSE,TemperatureRMSE_C, ...
    HumidityRatioRMSE,CO2RMSE_ppm,SafetyViolationSamples, ...
    TemperatureViolationSamples,HumidityViolationSamples, ...
    CO2ViolationSamples,DehumidificationSaturationAtChecks, ...
    AllDeviceSaturationAtChecks,InterlockConflicts, ...
    MaximumHumidityRatio,MaximumHumidityExcess,DeviceEnergy_kWh, ...
    ZeroSafetyViolation,NoDehumidificationSaturation, ...
    NoInterlockConflict,SafetyFeasible,CapacityFeasible);

feasibleIndex = find(CapacityFeasible,1,'first');
if isempty(feasibleIndex)
    selectedIndex = NaN;
    selectedUnitCount = NaN;
    selectedCapacity_kg_h = NaN;
    warning(['3、4、5台配置中没有同时实现零越界和零除湿饱和。', ...
        '请扩大容量候选，但不要直接修改安全边界。']);
else
    selectedIndex = feasibleIndex;
    selectedUnitCount = UnitCount(selectedIndex);
    selectedCapacity_kg_h = TotalCapacity_kg_h(selectedIndex);
end

capacityTable.Recommended = false(height(capacityTable),1);
if isfinite(selectedIndex)
    capacityTable.Recommended(selectedIndex) = true;
end

%% 8. 建立逐仓越界统计表
nRow = nCandidate*clock.nRoom;
CandidateUnitCount = repelem(UnitCount,clock.nRoom);
Room = repmat(cfg.room.names,nCandidate,1);
Stage = repmat(cfg.stage.names(cfg.room.stageIndex),nCandidate,1);
AllViolationSamples = reshape(roomViolationCount',nRow,1);
HumidityViolationSamplesByRoom = ...
    reshape(roomHumidityViolationCount',nRow,1);

roomViolationTable = table(CandidateUnitCount,Room,Stage, ...
    AllViolationSamples,HumidityViolationSamplesByRoom);

%% 9. 保存数据和所选配置的详细结果
writetable(capacityTable,fullfile(resultDirectory, ...
    'dehumidifier_capacity_scan.csv'));
writetable(roomViolationTable,fullfile(resultDirectory, ...
    'room_violation_summary.csv'));

selectedWeatherTable = table((0:71)',weather.Date,weather.T, ...
    weather.RH,weather.P,weather.w,weather.CO2ppm, ...
    'VariableNames',{'HourIndex','DateTime','T_out_C','RH_out_pct', ...
    'P_out_Pa','w_out_kg_kgda','CO2_out_ppm'});
writetable(selectedWeatherTable,fullfile(resultDirectory, ...
    'selected_weather_72h.csv'));

Parameter = {'unitCandidates';'perUnitCapacity'; ...
    'selectedWeatherStart';'selectionRule';'deviceVerification'};
Value = {mat2str(unitCandidates'); ...
    sprintf('%.15g',cfg.device.dehPerUnit_kg_h); ...
    char(selection.SelectedStart); ...
    ['选择零安全越界、零除湿饱和且容量裕量非负的', ...
     '最小每仓台数'];cfg.device.status};
Unit = {'台/仓';'kg/h';'北京时间';'-';'-'};
Category = {'算法设计值';'继承工况值';'派生选择值'; ...
    '算法设计规则';'待核对状态'};
writetable(table(Parameter,Value,Unit,Category), ...
    fullfile(resultDirectory,'parameter_provenance.csv'));

if isfinite(selectedIndex)
    cfgSelected = cfg;
    cfgSelected.device.nDehumidifier = selectedUnitCount;
    cfgSelected.limits.mWaterMin_kg_s = ...
        -(selectedCapacity_kg_h/3600);
    selectedCase = capacityCases{selectedIndex};
    selectedModel = model;
    selectedModel.device.nDehumidifier = selectedUnitCount;
    selectedModel.limits.mwMin = -(selectedCapacity_kg_h/3600);

    write_case_outputs(resultDirectory,'selected_capacity_periodic', ...
        selectedCase,clock,outdoor,cfgSelected,selectedModel);
    writetable(capacityTable(selectedIndex,:), ...
        fullfile(resultDirectory,'selected_capacity.csv'));
end

save(fullfile(resultDirectory,'dehumidifier_capacity_audit.mat'), ...
    'cfg','unitCandidates','capacityCases','capacityTable', ...
    'roomViolationTable','selectedIndex','selectedUnitCount', ...
    'selectedCapacity_kg_h','weather','selection','-v7.3');

%% 10. 绘制容量核查结果
fontName = choose_chinese_font();

fig1 = figure('Color','w','Position',[100 80 1050 650]);
yyaxis left;
b1 = bar(UnitCount,SafetyViolationSamples,0.35, ...
    'FaceColor',[0.10 0.45 0.75]);
ylabel('安全越界样本数');
yyaxis right;
b2 = bar(UnitCount+0.18,DehumidificationSaturationAtChecks,0.35, ...
    'FaceColor',[0.90 0.40 0.10]);
ylabel('除湿饱和次数');
xlabel('每仓除湿器台数');
title('不同除湿器配置下的安全越界与除湿饱和');
xticks(UnitCount); grid on; box on;
legend([b1 b2],{'安全越界样本数','除湿饱和次数'}, ...
    'Location','best');
set(findall(fig1,'-property','FontName'),'FontName',fontName);
exportgraphics(fig1,fullfile(resultDirectory, ...
    '图1_容量配置与越界饱和.png'),'Resolution',300);

fig2 = figure('Color','w','Position',[120 100 1000 620]);
plot(UnitCount,TotalCapacity_kg_h,'-o','LineWidth',1.8, ...
    'MarkerSize',7); hold on;
plot(UnitCount,PeakRequestedDehumidification_kg_h,'-s', ...
    'LineWidth',1.8,'MarkerSize',7);
yline(PeakRequestedDehumidification_kg_h(1), ...
    '--k','3台配置轨迹的最大需求派生值', ...
    'LineWidth',1.1);
xlabel('每仓除湿器台数');
ylabel('除湿量/(kg/h)');
title('设备总能力与控制器最大除湿需求');
xticks(UnitCount); grid on; box on;
legend({'设备总能力','本配置下最大请求量','原轨迹最大需求'}, ...
    'Location','best');
set(findall(fig2,'-property','FontName'),'FontName',fontName);
exportgraphics(fig2,fullfile(resultDirectory, ...
    '图2_除湿能力与最大需求.png'),'Resolution',300);

fig3 = figure('Color','w','Position',[140 60 1100 820]);
timeHour = clock.t/3600;
for q = 1:nCandidate
    subplot(nCandidate,1,q);
    plot(timeHour,capacityCases{q}.w(:,4),'LineWidth',1.1); hold on;
    plot(timeHour,capacityCases{q}.w(:,5),'LineWidth',1.1);
    yline(model.wMax(3),'--r','G3含湿量上限','LineWidth',1.1);
    grid on;
    ylabel('kg/kg干空气');
    title(sprintf('每仓%d台除湿器：4号仓与5号仓含湿量', ...
        UnitCount(q)));
    if q==nCandidate
        xlabel('时间/h');
    end
    legend({'4号仓','5号仓','G3安全上限'},'Location','best');
end
set(findall(fig3,'-property','FontName'),'FontName',fontName);
exportgraphics(fig3,fullfile(resultDirectory, ...
    '图3_G3仓含湿量容量对比.png'),'Resolution',300);

%% 11. 命令行结论
fprintf('\n====================================================\n');
fprintf('除湿设备容量可行性核查完成\n');
fprintf('====================================================\n');
disp(capacityTable);

if isfinite(selectedIndex)
    fprintf('\n最小可行配置：每仓%d台除湿器\n',selectedUnitCount);
    fprintf('总除湿能力 = %.6f kg/h\n',selectedCapacity_kg_h);
    fprintf('安全越界样本数 = %d\n', ...
        SafetyViolationSamples(selectedIndex));
    fprintf('除湿饱和次数 = %d\n', ...
        DehumidificationSaturationAtChecks(selectedIndex));
    fprintf('设备能耗 = %.6f kWh\n', ...
        DeviceEnergy_kWh(selectedIndex));
    fprintf('\n容量仿真通过后，可暂时将唯一参数入口改为：\n');
    fprintf('cfg.device.nDehumidifier = %d; %% [A] 容量扫描值\n', ...
        selectedUnitCount);
    fprintf(['注意：只有核对单台设备规格后，', ...
        '才能将cfg.device.verified设为true。\n']);
else
    fprintf('\n当前候选范围没有得到可行配置。\n');
    fprintf('应扩大台数或更换更大除湿能力的设备后重新核查。\n');
end

fprintf('\n结果目录：%s\n',resultDirectory);
fprintf('====================================================\n');
end

function model = build_model_from_cfg(cfg,H)
% 只进行字段映射，不重新手填物理参数。
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

function result = simulate_case(mode,thetaF,initial,model,clock,outdoor, ...
    eventNumericalTolerance)
% 周期组和固定触发组共用同一仿真函数，唯一差异为通信更新规则。
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
satQ = false(nStep,nRoom);
satW = false(nStep,nRoom);
satF = false(nStep,nRoom);
satDetail = false(nStep,nRoom,6);

isEventCheck = false(nStep,1);
isBroadcast = false(nStep,nRoom);
isFullControlExecution = false(nStep,nRoom);
eventErrorNorm = NaN(nStep,nRoom);
currentStateNorm = NaN(nStep,nRoom);
eventThreshold = NaN(nStep,nRoom);

zHat = zeros(nRoom,3);
hasInitialBroadcast = false;
uRawHold = zeros(nRoom,3);
uHold = zeros(nRoom,3);
satHold = false(nRoom,6);

for k = 1:nStep
    out = struct('T',outdoor.T(k),'w',outdoor.w(k), ...
        'P',outdoor.P(k),'Cppm',outdoor.Cppm(k), ...
        'c',outdoor.c(k));

    if mod(k-1,clock.nPhysicsPerCheck)==0
        isEventCheck(k) = true;
        zCurrent = normalized_state(T(k,:),w(k,:),c(k,:),model);

        if ~hasInitialBroadcast
            % 两组在t=0均强制发送，建立一致的初始网络状态。
            trigger = true(nRoom,1);
            errorNorm = zeros(nRoom,1);
            stateNorm = sqrt(sum(zCurrent.^2,2));
            thresholdNow = zeros(nRoom,1);
            hasInitialBroadcast = true;
        elseif strcmpi(mode,'Periodic')
            trigger = true(nRoom,1);
            errorNorm = sqrt(sum((zHat-zCurrent).^2,2));
            stateNorm = sqrt(sum(zCurrent.^2,2));
            thresholdNow = zeros(nRoom,1);
        elseif strcmpi(mode,'Fixed')
            eventError = zHat-zCurrent;
            errorNorm = sqrt(sum(eventError.^2,2));
            stateNorm = sqrt(sum(zCurrent.^2,2));
            thresholdNow = thetaF*stateNorm;

            % 第一项是正式固定触发条件；第二项只排除浮点0>=0伪触发。
            trigger = errorNorm>=thresholdNow & ...
                errorNorm>eventNumericalTolerance;
        else
            error('未知仿真模式：%s',mode);
        end

        zHat(trigger,:) = zCurrent(trigger,:);
        isBroadcast(k,:) = trigger';
        eventErrorNorm(k,:) = errorNorm';
        currentStateNorm(k,:) = stateNorm';
        eventThreshold(k,:) = thresholdNow';

        % 某仓发送后，它自己和两个邻居的协同反馈项需要刷新。
        affected = trigger | (model.A*double(trigger)>0);
        isFullControlExecution(k,:) = affected';

        % =============================================================
        % 固定事件触发控制律第一次实际使用
        % =============================================================
        % zHat保存各仓最近一次发送的归一化状态偏差。
        % etaHat = zHat + gamma*L*zHat：
        % 第一项纠正本仓相对阶段目标的偏差，第二项抑制相邻仓归一化
        % 偏差不一致。随后用H_g把etaHat转换为热量、水分和通风修正，
        % 再与由当前逐时室外状态计算的平衡前馈输入相加。
        % 未触发的仓不会更新zHat，因此反馈部分继续使用最近发送状态。
        [uRawHold,uHold,satHold] = ...
            controller_from_transmitted_state(zHat,model,out);
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

result = struct;
result.Name = mode;
result.thetaF = thetaF;
result.T = T;
result.w = w;
result.c = c;
result.Qraw = Qraw;
result.mwRaw = mwRaw;
result.Fraw = Fraw;
result.Qcmd = Qcmd;
result.mwCmd = mwCmd;
result.Fcmd = Fcmd;
result.satQ = satQ;
result.satW = satW;
result.satF = satF;
result.satDetail = satDetail;
result.isEventCheck = isEventCheck;
result.isBroadcast = isBroadcast;
result.isFullControlExecution = isFullControlExecution;
result.eventErrorNorm = eventErrorNorm;
result.currentStateNorm = currentStateNorm;
result.eventThreshold = eventThreshold;
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
    controller_from_transmitted_state(zHat,model,out)
etaHat = zHat+model.gamma*(model.L*zHat);
uFF = equilibrium_feedforward(model,out);
nRoom = numel(model.roomStage);
uRaw = zeros(nRoom,3);
for i = 1:nRoom
    g = model.roomStage(i);
    uRaw(i,:) = uFF(g,:)+(model.H(:,:,g)*etaHat(i,:)')';
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

function result = analyze_case(result,model,clock)
nStep = clock.nStep;
nRoom = clock.nRoom;
Cppm = mass_to_ppm(result.c,result.T, ...
    model.par.PIndoorRef_Pa,model.par);

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
    targetT(i) = model.T0(g);
    targetW(i) = model.w0(g);
    targetC(i) = model.Cppm0(g);
    safeTMin(i) = model.TMin(g);
    safeTMax(i) = model.TMax(g);
    safeWMin(i) = model.wMin(g);
    safeWMax(i) = model.wMax(g);
    safeCMin(i) = model.CppmMin(g);
    safeCMax(i) = model.CppmMax(g);
    normT(i) = model.scaleT(g);
    normW(i) = model.scaleW(g);
    normC(i) = model.scaleCppm(g);
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

[cool,heat,humid,dehumid,~,~,~,~,~,~,~,devicePower] = ...
    decode_devices(result.Qcmd,result.mwCmd,result.Fcmd,model.device);

interlockConflicts = nnz(cool>1e-12 & heat>1e-12) ...
    +nnz(humid>1e-12 & dehumid>1e-12);

checkIndex = find(result.isEventCheck);
nBroadcast = nnz(result.isBroadcast);
nReceive = sum(model.A,'all')/nRoom*nBroadcast;
activeCommunicationEnergy_J = ...
    (0.366*nBroadcast+0.310*nReceive)/1000;
deviceEnergy_kWh = sum(devicePower(1:end-1,:),'all') ...
    *clock.dt/3600;

metric = struct;
metric.NormalizedTrackingRMSE = sqrt(mean(z.^2,'all'));
metric.TemperatureRMSE_C = sqrt(mean(dT.^2,'all'));
metric.HumidityRatioRMSE = sqrt(mean(dW.^2,'all'));
metric.CO2RMSE_ppm = sqrt(mean(dC.^2,'all'));
metric.MaxAbsNormalizedError = max(abs(z),[],'all');
metric.SafetyViolationSamples = nnz(anyViolation);
metric.TemperatureViolationSamples = nnz(violT);
metric.HumidityViolationSamples = nnz(violW);
metric.CO2ViolationSamples = nnz(violC);
metric.SaturationAtChecks = nnz(result.satQ(checkIndex,:)) ...
    +nnz(result.satW(checkIndex,:)) ...
    +nnz(result.satF(checkIndex,:));
metric.InterlockConflicts = interlockConflicts;
metric.EventChecks = numel(checkIndex)*nRoom;
metric.NetworkTransmissions = nBroadcast;
metric.NetworkReceptions = nReceive;
metric.FullControlExecutions = nnz(result.isFullControlExecution);
metric.ActiveCommunicationEnergy_J = activeCommunicationEnergy_J;
metric.DeviceEnergy_kWh = deviceEnergy_kWh;

result.Cppm = Cppm;
result.dT = dT;
result.dW = dW;
result.dC = dC;
result.z = z;
result.anyViolation = anyViolation;
result.violT = violT;
result.violW = violW;
result.violC = violC;
result.devicePower_kW = devicePower;
result.metric = metric;
end

function scanTable = build_scan_table(candidateFraction,thetaCandidates, ...
    thetaStab,periodic,fixedCases,trackingTolerance_pct)
n = numel(thetaCandidates);
ThetaFraction = candidateFraction;
ThetaF = thetaCandidates;
ThetaStab = repmat(thetaStab,n,1);
NormalizedTrackingRMSE = zeros(n,1);
TrackingDegradation_pct = zeros(n,1);
TemperatureRMSE_C = zeros(n,1);
HumidityRatioRMSE = zeros(n,1);
CO2RMSE_ppm = zeros(n,1);
NetworkTransmissions = zeros(n,1);
NetworkReceptions = zeros(n,1);
CommunicationSaving_pct = zeros(n,1);
FullControlExecutions = zeros(n,1);
ControlExecutionSaving_pct = zeros(n,1);
SafetyViolationSamples = zeros(n,1);
TemperatureViolationSamples = zeros(n,1);
HumidityViolationSamples = zeros(n,1);
CO2ViolationSamples = zeros(n,1);
SaturationAtChecks = zeros(n,1);
InterlockConflicts = zeros(n,1);
ActiveCommunicationEnergy_J = zeros(n,1);
DeviceEnergy_kWh = zeros(n,1);

for q = 1:n
    m = fixedCases{q}.metric;
    NormalizedTrackingRMSE(q) = m.NormalizedTrackingRMSE;
    TrackingDegradation_pct(q) = 100* ...
        (m.NormalizedTrackingRMSE- ...
         periodic.metric.NormalizedTrackingRMSE) ...
        /periodic.metric.NormalizedTrackingRMSE;
    TemperatureRMSE_C(q) = m.TemperatureRMSE_C;
    HumidityRatioRMSE(q) = m.HumidityRatioRMSE;
    CO2RMSE_ppm(q) = m.CO2RMSE_ppm;
    NetworkTransmissions(q) = m.NetworkTransmissions;
    NetworkReceptions(q) = m.NetworkReceptions;
    CommunicationSaving_pct(q) = 100*(1- ...
        m.NetworkTransmissions/periodic.metric.NetworkTransmissions);
    FullControlExecutions(q) = m.FullControlExecutions;
    ControlExecutionSaving_pct(q) = 100*(1- ...
        m.FullControlExecutions/periodic.metric.FullControlExecutions);
    SafetyViolationSamples(q) = m.SafetyViolationSamples;
    TemperatureViolationSamples(q) = m.TemperatureViolationSamples;
    HumidityViolationSamples(q) = m.HumidityViolationSamples;
    CO2ViolationSamples(q) = m.CO2ViolationSamples;
    SaturationAtChecks(q) = m.SaturationAtChecks;
    InterlockConflicts(q) = m.InterlockConflicts;
    ActiveCommunicationEnergy_J(q) = m.ActiveCommunicationEnergy_J;
    DeviceEnergy_kWh(q) = m.DeviceEnergy_kWh;
end

WithinTrackingTolerance = ...
    TrackingDegradation_pct<=trackingTolerance_pct;

% 绝对安全：固定触发72 h内一次安全越界也没有。
AbsoluteSafetyPass = SafetyViolationSamples==0;

% 相对安全：固定触发没有比周期通信基准增加安全越界。
% 该指标只用于隔离事件触发机制的增量影响，不能代替绝对安全。
NoAdditionalSafetyViolation = SafetyViolationSamples ...
    <= periodic.metric.SafetyViolationSamples;

InsideStabilityBound = ThetaF<ThetaStab;
NoInterlockConflict = InterlockConflicts==0;

QualifiedForCalibration = WithinTrackingTolerance ...
    & NoAdditionalSafetyViolation ...
    & InsideStabilityBound & NoInterlockConflict;

QualifiedForFinalSafeUse = QualifiedForCalibration ...
    & AbsoluteSafetyPass;

% 保留旧字段名，兼容可能调用该字段的后处理代码。
Qualified = QualifiedForCalibration;

scanTable = table(ThetaFraction,ThetaF,ThetaStab, ...
    NormalizedTrackingRMSE,TrackingDegradation_pct, ...
    TemperatureRMSE_C,HumidityRatioRMSE,CO2RMSE_ppm, ...
    NetworkTransmissions,NetworkReceptions,CommunicationSaving_pct, ...
    FullControlExecutions,ControlExecutionSaving_pct, ...
    SafetyViolationSamples,TemperatureViolationSamples, ...
    HumidityViolationSamples,CO2ViolationSamples, ...
    SaturationAtChecks,InterlockConflicts, ...
    ActiveCommunicationEnergy_J,DeviceEnergy_kWh, ...
    WithinTrackingTolerance,AbsoluteSafetyPass, ...
    NoAdditionalSafetyViolation,InsideStabilityBound, ...
    NoInterlockConflict,QualifiedForCalibration, ...
    QualifiedForFinalSafeUse,Qualified);
end

function selectedIndex = choose_best_threshold(scanTable,candidateIndex)
% 在指定候选集合中执行统一的三级择优规则：
% 1. 通信节省率最大；
% 2. 若并列，跟踪误差增幅最小；
% 3. 若仍并列，阈值较小。

assert(~isempty(candidateIndex), ...
    'choose_best_threshold收到的候选集合为空。');

saving = scanTable.CommunicationSaving_pct(candidateIndex);
maximumSaving = max(saving);
tie = candidateIndex(abs(saving-maximumSaving)<1e-12);

degradation = scanTable.TrackingDegradation_pct(tie);
minimumDegradation = min(degradation);
tie2 = tie(abs(degradation-minimumDegradation)<1e-12);

[~,localIndex] = min(scanTable.ThetaF(tie2));
selectedIndex = tie2(localIndex);
end

function tableOut = metric_to_table(result,condition)
m = result.metric;
Condition = {condition};
NormalizedTrackingRMSE = m.NormalizedTrackingRMSE;
TemperatureRMSE_C = m.TemperatureRMSE_C;
HumidityRatioRMSE = m.HumidityRatioRMSE;
CO2RMSE_ppm = m.CO2RMSE_ppm;
NetworkTransmissions = m.NetworkTransmissions;
NetworkReceptions = m.NetworkReceptions;
FullControlExecutions = m.FullControlExecutions;
SafetyViolationSamples = m.SafetyViolationSamples;
TemperatureViolationSamples = m.TemperatureViolationSamples;
HumidityViolationSamples = m.HumidityViolationSamples;
CO2ViolationSamples = m.CO2ViolationSamples;
SaturationAtChecks = m.SaturationAtChecks;
InterlockConflicts = m.InterlockConflicts;
ActiveCommunicationEnergy_J = m.ActiveCommunicationEnergy_J;
DeviceEnergy_kWh = m.DeviceEnergy_kWh;
tableOut = table(Condition,NormalizedTrackingRMSE,TemperatureRMSE_C, ...
    HumidityRatioRMSE,CO2RMSE_ppm,NetworkTransmissions, ...
    NetworkReceptions,FullControlExecutions,SafetyViolationSamples, ...
    TemperatureViolationSamples,HumidityViolationSamples, ...
    CO2ViolationSamples, ...
    SaturationAtChecks,InterlockConflicts, ...
    ActiveCommunicationEnergy_J,DeviceEnergy_kWh);
end

function write_case_outputs(resultDirectory,prefix,result,clock,outdoor, ...
    cfg,model)
nStep = clock.nStep;
nRoom = clock.nRoom;
timeLong = repmat(clock.t,nRoom,1);
roomIndexLong = kron((1:nRoom)',ones(nStep,1));
roomLong = cfg.room.names(roomIndexLong);
stageLong = cfg.stage.names(cfg.room.stageIndex(roomIndexLong));

stateTable = table(timeLong,timeLong/3600,roomLong,stageLong, ...
    result.T(:),result.w(:),result.Cppm(:), ...
    result.dT(:),result.dW(:),result.dC(:), ...
    reshape(result.z(:,:,1),[],1), ...
    reshape(result.z(:,:,2),[],1), ...
    reshape(result.z(:,:,3),[],1), ...
    result.violT(:),result.violW(:),result.violC(:), ...
    result.anyViolation(:), ...
    repmat(result.isEventCheck,nRoom,1), ...
    result.isBroadcast(:), ...
    result.eventErrorNorm(:),result.currentStateNorm(:), ...
    result.eventThreshold(:), ...
    repmat(outdoor.T,nRoom,1),repmat(outdoor.w,nRoom,1), ...
    repmat(outdoor.Cppm,nRoom,1), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','T_C', ...
    'w_kg_kgda','CO2_ppm','dT_C','dw_kg_kgda','dCO2_ppm', ...
    'zT','zw','zCO2','TemperatureViolation', ...
    'HumidityViolation','CO2Violation','AnySafetyViolation','EventCheck', ...
    'Broadcast','EventErrorNorm','CurrentStateNorm','EventThreshold', ...
    'T_out_C','w_out_kg_kgda','CO2_out_ppm'});
writetable(stateTable,fullfile(resultDirectory, ...
    [prefix '_state_trajectories.csv']));

[cool,heat,humid,dehumid,fanLow,fanHigh,fanOff,coolDuty, ...
    heatDuty,humidDuty,dehumidDuty,power] = decode_devices( ...
    result.Qcmd,result.mwCmd,result.Fcmd,model.device);

deviceTable = table(timeLong,timeLong/3600,roomLong,stageLong, ...
    result.Qraw(:),result.Qcmd(:),cool(:),heat(:), ...
    coolDuty(:),heatDuty(:),result.mwRaw(:),result.mwCmd(:), ...
    humid(:),dehumid(:),humidDuty(:),dehumidDuty(:), ...
    result.Fraw(:),result.Fcmd(:),result.Fcmd(:)*3600, ...
    fanOff(:),fanLow(:),fanHigh(:),power(:), ...
    result.satQ(:),result.satW(:),result.satF(:), ...
    'VariableNames',{'Time_s','Time_h','Room','Stage','Q_raw_W', ...
    'Q_cmd_W','Cooling_kW','Heating_kW','CoolingDuty','HeatingDuty', ...
    'mw_raw_kg_s','mw_cmd_kg_s','Humidification_kg_h', ...
    'Dehumidification_kg_h','HumidifierDuty','DehumidifierDuty', ...
    'F_raw_m3_s','F_cmd_m3_s','F_cmd_m3_h','FanOffDuty', ...
    'FanLowDuty','FanHighDuty','DevicePower_kW', ...
    'ThermalSaturation','MoistureSaturation','VentSaturation'});
writetable(deviceTable,fullfile(resultDirectory, ...
    [prefix '_device_commands.csv']));

violationTable = stateTable(result.anyViolation(:),:);
writetable(violationTable,fullfile(resultDirectory, ...
    [prefix '_constraint_violations.csv']));
end

function write_provenance(filePath,cfg,weatherFile,selection,thetaStab, ...
    candidateFraction,trackingTolerance,eventTolerance)
Parameter = {'weatherFile';'selectedWeatherStart'; ...
    'weatherSelectionRule';'dt';'eventCheckPeriod';'duration'; ...
    'gamma';'topology';'roomStageAssignment'; ...
    'room2InitialHumidityRatio';'thetaStab'; ...
    'candidateFraction';'trackingTolerance'; ...
    'eventNumericalTolerance';'deviceVerificationStatus'};
Value = {weatherFile;char(selection.SelectedStart); ...
    ['72 h daily windows; annual-mean target for means; ', ...
     'candidate-median target for variations; IQR scaling']; ...
    sprintf('%g',cfg.time.dt_s); ...
    sprintf('%g',cfg.time.eventCheckPeriod_s); ...
    sprintf('%g',cfg.time.mainDuration_s); ...
    sprintf('%.15g',cfg.control.gamma); ...
    'undirected ring 1-2-3-4-5-6-1';'[G1,G1,G2,G3,G3,G4]'; ...
    '0.0165089';sprintf('%.15g',thetaStab); ...
    mat2str(candidateFraction',4); ...
    sprintf('%.15g',trackingTolerance); ...
    sprintf('%.15g',eventTolerance);cfg.device.status};
Unit = {'path';'datetime';'-';'s';'s';'s';'-';'-';'-'; ...
    'kg/kgda';'-';'fraction';'%';'-';'-'};
Category = {'用户提供值';'派生计算值';'算法设计值'; ...
    '算法设计值';'算法设计值';'算法设计值'; ...
    '算法设计值';'算法设计值';'用户指定值';'用户指定值'; ...
    '派生计算值';'算法设计值';'算法设计值'; ...
    '算法实现容差';'待核对值'};
SourceNote = {'用户提供的南宁2023逐时重建气象表'; ...
    '按固定规则自动选择，非按控制结果挑选'; ...
    '与实验二保持一致';'唯一参数入口';'唯一参数入口'; ...
    '唯一参数入口';'实验二敏感性分析确定'; ...
    '六仓协同控制拓扑';'用户确定的阶段分配'; ...
    '明确排除作废值0.0171547'; ...
    '全精度Lyapunov复算结果'; ...
    '理论上界5%至90%，间隔5%'; ...
    '本研究控制性能保持条件，不是文献标准'; ...
    '仅防止浮点0>=0伪触发'; ...
    '继承当前设备参数，仍需完成厂商资料核对'};
writetable(table(Parameter,Value,Unit,Category,SourceNote),filePath);
end

function make_figures(resultDirectory,scanTable,periodic,recommended, ...
    clock,cfg,thetaStab,trackingTolerance)
fontName = choose_chinese_font();
visible = 'on';

fig1 = figure('Visible',visible,'Color','w', ...
    'Position',[100 80 1080 720]);
yyaxis left;
plot(scanTable.ThetaF,scanTable.TrackingDegradation_pct, ...
    '-o','LineWidth',1.4,'MarkerSize',5);
hold on;
yline(trackingTolerance,'--r','1%控制性能容差','LineWidth',1.2);
ylabel('综合跟踪误差增幅/%');
yyaxis right;
plot(scanTable.ThetaF,scanTable.CommunicationSaving_pct, ...
    '-s','LineWidth',1.4,'MarkerSize',5);
ylabel('通信节省率/%');
xline(thetaStab,':k','理论稳定上界','LineWidth',1.2);
xlabel('固定事件触发阈值 \theta_f');
title('固定阈值对控制性能与通信节省的影响');
grid on;
set(findall(fig1,'-property','FontName'),'FontName',fontName);
exportgraphics(fig1,fullfile(resultDirectory, ...
    '图1_固定阈值扫描主结果.png'),'Resolution',180);

fig2 = figure('Visible',visible,'Color','w', ...
    'Position',[120 100 1080 620]);
bar(scanTable.ThetaF,scanTable.NetworkTransmissions);
hold on;
yline(periodic.metric.NetworkTransmissions,'--r', ...
    '周期通信基准','LineWidth',1.2);
xlabel('固定事件触发阈值 \theta_f');
ylabel('全网发送次数');
title('不同固定阈值下的通信次数');
grid on;
set(findall(fig2,'-property','FontName'),'FontName',fontName);
exportgraphics(fig2,fullfile(resultDirectory, ...
    '图2_不同阈值通信次数.png'),'Resolution',180);

if isempty(recommended)
    return;
end

timeHour = clock.t/3600;
fig3 = figure('Visible',visible,'Color','w', ...
    'Position',[140 80 1160 800]);
subplot(3,2,1);
plot(timeHour,periodic.T,'LineWidth',0.8); grid on;
ylabel('温度/℃'); title('周期通信：温度');
subplot(3,2,2);
plot(timeHour,recommended.T,'LineWidth',0.8); grid on;
ylabel('温度/℃'); title('初步标定固定阈值：温度');
subplot(3,2,3);
plot(timeHour,periodic.w,'LineWidth',0.8); grid on;
ylabel('含湿量/(kg/kg干空气)'); title('周期通信：含湿量');
subplot(3,2,4);
plot(timeHour,recommended.w,'LineWidth',0.8); grid on;
ylabel('含湿量/(kg/kg干空气)'); title('初步标定固定阈值：含湿量');
subplot(3,2,5);
plot(timeHour,periodic.Cppm,'LineWidth',0.8); grid on;
ylabel('CO₂/ppm'); xlabel('时间/h'); title('周期通信：CO₂');
subplot(3,2,6);
plot(timeHour,recommended.Cppm,'LineWidth',0.8); grid on;
ylabel('CO₂/ppm'); xlabel('时间/h'); title('初步标定固定阈值：CO₂');
sgtitle(sprintf('周期通信与初步标定固定阈值状态轨迹对比（thetaF=%.6g）', ...
    recommended.thetaF));
set(findall(fig3,'-property','FontName'),'FontName',fontName);
exportgraphics(fig3,fullfile(resultDirectory, ...
    '图3_周期与初步固定阈值状态对比.png'),'Resolution',180);

fig4 = figure('Visible',visible,'Color','w', ...
    'Position',[160 120 1120 620]);
checkIndex = find(recommended.isEventCheck);
imagesc(timeHour(checkIndex),1:clock.nRoom, ...
    double(recommended.isBroadcast(checkIndex,:)'));
colormap([1 1 1;0.05 0.35 0.75]);
caxis([0 1]);
colorbar('Ticks',[0 1],'TickLabels',{'未触发','触发'});
xlabel('时间/h'); ylabel('仓号');
yticks(1:clock.nRoom); yticklabels(cfg.room.names);
title(sprintf('初步标定固定阈值的六仓事件触发时刻（thetaF=%.6g）', ...
    recommended.thetaF));
set(findall(fig4,'-property','FontName'),'FontName',fontName);
exportgraphics(fig4,fullfile(resultDirectory, ...
    '图4_初步固定阈值事件触发分布.png'),'Resolution',180);
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
[bestScore,bestCandidate] = min(score);
idx = startIndex(bestCandidate):(startIndex(bestCandidate)+71);

weather = struct('Date',date(idx),'T',T(idx),'RH',RH(idx), ...
    'P',P(idx),'CO2ppm',CO2ppm(idx),'w',w(idx));
selection = struct('SelectedStart',date(idx(1)), ...
    'SelectedEnd',date(idx(end)), ...
    'CandidateCount',nCandidate,'Score',bestScore);
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
