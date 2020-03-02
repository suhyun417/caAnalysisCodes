function [] = saveTimingML_DFL(fname_bhv, fname_out_mat)
% get behavioral data and timing info from ML

[data, MLConfig, TrialRecord] = mlread(fname_bhv);
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

t.trialStart = data.AbsoluteTrialStartTime;
t.sendTTL_start = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==900);
t.startCa = find(data.AnalogData.Button.Btn1>0, 1);
t.reward = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==90);
t.sendTTL_end = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==990);
t.movieOnset = data.BehavioralCodes.CodeTimes(data.BehavioralCodes.CodeNumbers==40);

% catCodeTimes = cat(1, data.BehavioralCodes.CodeTimes);
% catCodeNumbers = cat(1, data.BehavioralCodes.CodeNumbers);

t_adj.trialStart = t.trialStart - t.startCa;
t_adj.sendTTL_start = t.sendTTL_start - t.startCa;
t_adj.sendTTL_end = t.sendTTL_end - t.startCa;
t_adj.movieOnset = t.movieOnset - t.startCa;
t_adj.reward = t.reward - t.startCa;

analog.sampleInterval = data.AnalogData.SampleInterval;
analog.photodiode = data.AnalogData.PhotoDiode;
analog.eye = cat(2, data.AnalogData.Eye, data.AnalogData.General.Gen1);
analog.CaSync = data.AnalogData.Button.Btn1;

save(fname_out_mat, 't', 't_adj', 'stim', 'analog', 'MLConfig', 'TrialRecord')