% temp_analyzeML.m
%
% analyze bhv2 files from Monkeylogic for BPM
% quick and dirty script to try out stuff

dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %'/Users/parks20/ResearchProjects/0MARMOSET/0FaceNeurons/_data/invivoCalciumImaging/Max/20191010_Max'; % 20191004_Max/';

nameSubj = 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
dateSession = '20191113'; %'20191125';  

dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);


% 191125_Max_Ca_BPM_122309.bhv2
d_BPM = dir(fullfile(dirBHV, sprintf('%s_%s_Ca_BPM_*.bhv2', datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') , nameSubj)));
d_DFL = dir(fullfile(dirBHV, sprintf('%s_%s_Ca_DFL_*.bhv2', datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') , nameSubj)));

d_BPM = d_BPM(2:3);
nRun = length(d_BPM);

for iRun = 1:nRun

    saveFileName = fullfile(dirProcdata_session, sprintf('BPM_%d_tML.mat', iRun));
    DataFile = fullfile(dirBHV, d_BPM(iRun).name);
    
    [data_cat,MLConfig,TrialRecord] = mlconcatenate_woUserVars(DataFile);
    t.sendTTL = data_cat.BehavioralCodes(1).CodeTimes(2);
    t.startCa = find(data_cat.AnalogData.Button.Btn1>0, 1);
    t.trialStart = data_cat.AbsoluteTrialStartTime;
    
    catCodeTimes = cat(1, data_cat.BehavioralCodes.CodeTimes);
    catCodeNumbers = cat(1, data_cat.BehavioralCodes.CodeNumbers);
    t.blankOnset_beforeStim = catCodeTimes(catCodeNumbers==15);
    t.stimOnset = catCodeTimes(catCodeNumbers==20);
    t.blankOnset_afterStim = catCodeTimes(catCodeNumbers==25);
    
    t_adj.trialStart = t.trialStart - t.startCa;
    t_adj.blankOnset_beforeStim = t.blankOnset_beforeStim - t.startCa;
    t_adj.stimOnset = t.stimOnset - t.startCa;
    t_adj.blankOnset_afterStim = t.blankOnset_afterStim - t.startCa;
    
    data = mlread(DataFile);
    stim.condMat = data(1).UserVars.UserDefined.condMat;
    stim.trialError = data_cat.TrialError;       
    
    t.sig.CaSync = data_cat.AnalogData.Button.Btn1;
    
    save(saveFileName, 't', 't_adj', 'stim', 'data_cat', 'MLConfig', 'TrialRecord', 'DataFile')

end


%% DFL
d_DFL = dir(fullfile(dirBHV, sprintf('%s_%s_Ca_DFL_*.bhv2', datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') , nameSubj)));

nRun = length(d_DFL);

for iRun = 1:nRun

    saveFileName = fullfile(dirProcdata_session, sprintf('DFL_%d_tML.mat', iRun));
    DataFile = fullfile(dirBHV, d_DFL(iRun).name);
    
    [data, MLConfig, TrialRecord] = mlread(DataFile);
    stim.cond = data.Condition;
    temp = data.TaskObject.Attribute{2};
    [p, name, ext] = fileparts(temp{2});
    stim.nameMovie = char(regexp(name, '(set\w{3})', 'match'));
    
    % Behavioral Codes
    % FP_ON = 20;
    % WAIT_FOR_TR = 30; % I made this tracker "null" so in theory it
    % doesn't wait for the trigger from camera 
    % MOVIE_ON = 40;
    % REWARD = 90;
    
%     [data_cat,MLConfig,TrialRecord] = mlconcatenate_woUserVars(DataFile);
    t.trialStart = data.AbsoluteTrialStartTime;
    t.sendTTL_start = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==900);
    t.startCa = find(data.AnalogData.Button.Btn1>0, 1);    
    t.reward = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==90);
    t.sendTTL_end = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==990);
    t.movieOnset = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==40);
    
    catCodeTimes = cat(1, data.BehavioralCodes.CodeTimes);
    catCodeNumbers = cat(1, data.BehavioralCodes.CodeNumbers);
%     t.blankOnset_beforeStim = catCodeTimes(catCodeNumbers==15);
%     t.stimOnset = catCodeTimes(catCodeNumbers==20);
%     t.blankOnset_afterStim = catCodeTimes(catCodeNumbers==25);
    
%     t_adj.trialStart = t.trialStart - t.startCa;
%     t_adj.sendTTL_start = t.sendTTL_start - t.startCa;  
    t_adj.sendTTL_end = t.sendTTL_end - t.startCa;
    t_adj.movieOnset = t.movieOnset - t.startCa;
    t_adj.reward = t.reward - t.startCa;
    
%     data = mlread(DataFile);
%     stim.condMat = data(1).UserVars.UserDefined.condMat;
%     stim.trialError = data_cat.TrialError;       
    
    t.sig.CaSync = data.AnalogData.Button.Btn1;
    
    save(saveFileName, 't', 't_adj', 'stim', 'data', 'MLConfig', 'TrialRecord', 'DataFile')

end

%% Read .bhv2 files (from MonkeyLogic) to get the timing information

% % for iF = 1:length(setFileName)
% iF = 2;
%     [data,MLConfig,TrialRecord] = mlconcatenate_woUserVars(fullfile(dirBHV, setFileName{iF}));
%     t_sendTTL = data.BehavioralCodes(1).CodeTimes(2);
%     t_startCa = find(data.AnalogData.Button.Btn1>0, 1);
%     dT = t_startCa-t_sendTTL;
%     fprintf(1, '         Filename: %s, t_sendTTL: %3.3f, t_startCa: %3.3f, dT: %3.3f\n', setFileName{iF}, t_sendTTL, t_startCa, dT);
% % end
% 
% 
% t_endTTL = data.BehavioralCodes(40).CodeTimes(data.BehavioralCodes(40).CodeNumbers==990);
% t_trialStart_abs = data.AbsoluteTrialStartTime;
% catCodeTimes = cat(1, data.BehavioralCodes.CodeTimes);
% catCodeNumbers = cat(1, data.BehavioralCodes.CodeNumbers);
% t_stimOnset = catCodeTimes(catCodeNumbers==20);
% t_blankOnset = catCodeTimes(catCodeNumbers==15);
% t_afterOff = catCodeTimes(catCodeNumbers==25)+3000;

% t_onset = t_stimOnset; %t_blankOnset;
% t_end = t_afterOff;

fs = 15; %10;
delay = 0; %0.2; %t_startCa/1000; %0.2;
t_onset_adj = t_onset/1000-delay;
locStart = floor(t_onset_adj./(1/fs));
locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));

clear matTS_trial
for iTrial = 1:length(locStart)
    matTS_trial{iTrial} = matTS(locStart(iTrial):locStart(iTrial)+49, :); %locEnd(iTrial), :); %locOnset_end(iTrial), :);
end

cond = cat(1, data.Condition(data.TrialError>1));
[idCond, locTrial] = sort(cond);

setCond = unique(idCond);
nCond = length(setCond);

% clear matTS_sorted %= [];
matTS_sorted = [];
for iCond = 1:nCond
    curLoc = locTrial(idCond==setCond(iCond));
    catCond = cat(3, matTS_trial{curLoc});
    matTS_sorted = cat(4, matTS_sorted, catCond); % time x cell x trials x conditions
end
