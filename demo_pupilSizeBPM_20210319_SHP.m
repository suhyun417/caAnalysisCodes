% demo_pupilSizeBPM_20210319_SHP.m
%
% To compute pupil size during stimulus presentation for each trial

%% Define filename. 
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have

%% Read the file
data = mlread(filename);

%% For each trial
iTrial = 1; % which trial?

% retrieve pupil data of this trial
pData = data(iTrial).AnalogData.General.Gen1; % pupil data

% get the time information from BehavioralCodes
% Note for CodeNumbers: 20 (stimulus on), 55 (stimulus off)
index_stimOn = find(data(iTrial).BehavioralCodes.CodeNumbers == 20);
index_stimOff = find(data(iTrial).BehavioralCodes.CodeNumbers == 55);

time_stimOn = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOn));
time_stimOff = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOff));

% stack the time information for each trial
results_time(iTrial, 1) = time_stimOn;
results_time(iTrial, 2) = time_stimOff;

% get pupil size during stimulus presentation (between stim on and stim off
% timing)
results_pupil(iTrial).pupilData_raw = pData(time_stimOn:time_stimOff);
results_pupil(iTrial).pupilData_mean = mean(pData(time_stimOn:time_stimOff));
results_pupil(iTrial).pupilData_ste = std(pData(time_stimOn:time_stimOff))./sqrt(numel(pData(time_stimOn:time_stimOff))-1);

%% To do the above computation for all the trials, we need to make it to a for-loop
%% "help for" will give you the basic syntax
%%  If the for loop was done correctly, the "results_time" and "results_pupil" will have 51 entries

