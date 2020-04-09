function [] = doSaveTimingML_BPM(fname_bhv, fname_out_mat)
% get behavioral data and timing info from ML 

%     saveFileName = fullfile(dirProcdata_session, sprintf('BPM_%d_tML.mat', iRun));
%     DataFile = fullfile(dirBHV, d_BPM(iRun).name);

[data_cat,MLConfig,TrialRecord] = mlconcatenate_woUserVars(fname_bhv);

infoTrial.trialError = data_cat.TrialError;
infoTrial.flagValidTrial = data_cat.TrialError>1;
infoTrial.indValidTrial = find(data_cat.TrialError>1);
infoTrial.numValidTrial = sum(infoTrial.flagValidTrial);

t.sendTTL = data_cat.BehavioralCodes(1).CodeTimes(2);
t.startCa = find(data_cat.AnalogData.Button.Btn1>0, 1);
t.trialStart = data_cat.AbsoluteTrialStartTime(infoTrial.flagValidTrial);

catCodeTimes = cat(1, data_cat.BehavioralCodes.CodeTimes);
catCodeNumbers = cat(1, data_cat.BehavioralCodes.CodeNumbers);
t.blankOnset_beforeStim = catCodeTimes(catCodeNumbers==15);
t.stimOnset = catCodeTimes(catCodeNumbers==20);
t.blankOnset_afterStim = catCodeTimes(catCodeNumbers==25);

t_adj.trialStart = t.trialStart - t.startCa;
t_adj.blankOnset_beforeStim = t.blankOnset_beforeStim - t.startCa;
t_adj.stimOnset = t.stimOnset - t.startCa;
t_adj.blankOnset_afterStim = t.blankOnset_afterStim - t.startCa;

analog.sampleInterval = data_cat.AnalogData.SampleInterval;
analog.photodiode = data_cat.AnalogData.PhotoDiode;
analog.eye = cat(2, data_cat.AnalogData.Eye, data_cat.AnalogData.General.Gen1);
analog.CaSync = data_cat.AnalogData.Button.Btn1;

% stimulus-related
data = mlread(fname_bhv);
infoTrial.condMat = data(infoTrial.indValidTrial(1)).UserVars.UserDefined.condMat(infoTrial.indValidTrial,:);
infoTrial.idStim = infoTrial.condMat(:,1).*10 + infoTrial.condMat(:,2); % unique ID for each stimulus item
[setIdStim, ia, ~] = unique(infoTrial.idStim);

condName = {'human face', 'marmoset face', 	'marmoset body', 'scene', 'non familiar object', 'familiar object',...
    'phase scrambled', 'space scrambled', 'grating', 'random dot motion'};

for iStim = 1:length(setIdStim)
    curStim = setIdStim(iStim);

    infoStim(iStim).idStim = curStim;
    infoStim(iStim).nameCondition = condName{floor(curStim/10)};
    if floor(curStim/10) < 9
        infoStim(iStim).stimulus = data_cat.ObjectStatusRecord(infoTrial.indValidTrial(ia(iStim))).SceneParam(3).AdapterArgs{2}{4,2}{1};
    else
        infoStim(iStim).stimulus = data_cat.ObjectStatusRecord(infoTrial.indValidTrial(ia(iStim))).SceneParam(3).AdapterArgs{2};
    end
end

infoTrial.uniqueSetIdStim = setIdStim;
infoTrial.infoStim = infoStim;

save(fname_out_mat, 't', 't_adj', 'infoTrial', 'analog', 'MLConfig', 'TrialRecord')

