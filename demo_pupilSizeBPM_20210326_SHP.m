% demo_pupilSizeBPM_20210326_SHP.m
%
% To compute pupil size during stimulus presentation for each trial

%% Define filename. 
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have

%% Read the file
data = mlread(filename);

%% For each trial
results_time = NaN(length(data), 2);
for iTrial = 1:length(data) % which trial?

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

end

% concatenate the pupil size summary statistics from each trial across trials
%: excluding the first trial, which is repeated as the last trial during the session
catPupilMean = cat(1, results_pupil(2:51).pupilData_mean); % average (mean) pupil size
catPupilSte = cat(1, results_pupil(2:51).pupilData_ste); % standard error of the mean

% quick histogram
figure;
histogram(catPupilMean)
histogram(catPupilMean, 30)


%% Condition information
trial_conditions = cat(1, data(2:51).Condition); 

% condition info: sometimes you need to manually hard-code things like this..
% semi-colon inbetween items make this as a column vector (i.e.
% concatenation along the first dimension): see "help cat"
condName = {'human face'; 'marmoset face';	'marmoset body'; 'scene'; 'non familiar object'; 'familiar object';...
    'phase scrambled'; 'space scrambled'; 'grating'; 'random dot motion'};

% % assessing data in cell array
% condName{1} % this gives you the content "inside" of that cell
% condName(1) % this gives you the cell itself

%% Draw two histograms of average pupil size, one for marmoset face conditions, and another for non-familiar object



%% Which trials were "marmoset face" trials?
condition_id = 2; % for marmoset face (refer to the hard-coded "condName" above)
trial_id_marmosetFace = find(trial_conditions == condition_id);

% you can put number directly as below
trial_id_marmosetFace = find(trial_conditions == 2);




