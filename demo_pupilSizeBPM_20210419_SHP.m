% demo_pupilSizeBPM_20210416_SHP.m
%
% To compute pupil size during stimulus presentation for each trial
% 

%% Define filename. 
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
nameSubj =  'Max'; %'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
% dateSession = '20191113'; %'20191125';  
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);
S = table2struct(infoSession);

% setExpName = {S.ExpName}';
setMLFilename = {S.MLFilename}';

indBPMRuns = contains(setMLFilename, 'BPM') & cat(1, S.flagPreproc)>0; %% containing "BPM" in filename AND flagPreproc value of 1
setFilename = setMLFilename(indBPMRuns);

for iFile = 1:length(setFilename)
filename = strcat(setFilename{iFile}, '.bhv2');

dateSession = filename(1:6);

if str2num(dateSession) < 191121
    dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
else
    dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
end

% filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have

%% Read the file
data = mlread(fullfile(dirBHV, filename)); % mlread(filename);

%% For each trial
% checking valid trials
locValidTrials = find(cat(1, data.TrialError)>1); % excluding the first trial becuase it's repeated as the last trial
if sum(ismember(locValidTrials, 1))>0
    locValidTrials = locValidTrials(~ismember(locValidTrials,1));
end

clear results*
results_time = NaN(length(locValidTrials), 2);

for iTT = 1:length(locValidTrials) % which trial?
    
    iTrial = locValidTrials(iTT);
    
    % retrieve pupil data of this trial
    pData = data(iTrial).AnalogData.General.Gen1; % pupil data
    
    % get the time information from BehavioralCodes
    % Note for CodeNumbers: 20 (stimulus on), 55 (stimulus off)
    index_stimOn = find(data(iTrial).BehavioralCodes.CodeNumbers == 20);
    index_stimOff = find(data(iTrial).BehavioralCodes.CodeNumbers == 55);
    
    time_stimOn = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOn));
    time_stimOff = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOff));
    
    % stack the time information for each trial
    results_time(iTT, 1) = time_stimOn;
    results_time(iTT, 2) = time_stimOff;
    
    % get pupil size during stimulus presentation (between stim on and stim off
    % timing)
    results_pupil(iTT).pupilData_raw = pData(time_stimOn:time_stimOff);
    results_pupil(iTT).pupilData_mean = mean(pData(time_stimOn:time_stimOff));
    results_pupil(iTT).pupilData_ste = std(pData(time_stimOn:time_stimOff))./sqrt(numel(pData(time_stimOn:time_stimOff))-1);

end

% concatenate the pupil size summary statistics from each trial across trials
catPupilMean=[]; catPupilSte=[];
catPupilMean = cat(1, results_pupil.pupilData_mean); % average (mean) pupil size
catPupilSte = cat(1, results_pupil.pupilData_ste); % standard error of the mean

% % quick histogram
% figure;
% histogram(catPupilMean)
% histogram(catPupilMean, 30)

Results_trial.results_pupil = results_pupil;
Results_trial.results_time = results_time;
Results_trial.info = ['results_pupil contains raw pupil size data, average and standard error for each trial.', 'results_time contains stimulus on (1st column) and off (2nd column) time stamps for each trial']; 

%% Condition information
trial_conditions=[];
trial_conditions = cat(1, data(locValidTrials).Condition);

% condition info: sometimes you need to manually hard-code things like this..
% semi-colon inbetween items make this as a column vector (i.e.
% concatenation along the first dimension): see "help cat"
condName = {'human face'; 'marmoset face';	'marmoset body'; 'scene'; 'non familiar object'; 'familiar object';...
    'phase scrambled'; 'space scrambled'; 'grating'; 'random dot motion'};

% % assessing data in cell array
% condName{1} % this gives you the content "inside" of that cell
% condName(1) % this gives you the cell itself

%% Compare different conditions
conditionType = unique(trial_conditions); % extract the condition types used in this session

matPupilMean_Cond = cell([10, 1]); %{}; % [];
for condition_id = 1:10 %length(conditionType)
    if sum(ismember(conditionType, condition_id))>0
        trial_id = find(trial_conditions == condition_id);    
        matPupilMean_Cond{condition_id, 1} = catPupilMean(trial_id);
    else
        continue;
    end
%     matPupilMean_Cond(: , iCondition) = catPupilMean(trial_id);
end

Results_session(iFile).trial_conditions = trial_conditions;
Results_session(iFile).condName_all = condName;
Results_session(iFile).conditionType = conditionType;
Results_session(iFile).matPupilMean_Cond = matPupilMean_Cond;
Results_session(iFile).Results_trial = Results_trial;

% figure;
% plot(1:5, matPupilMean_Cond, 'bo')
% xlim([0 6])
% ylabel('Pupil Size')

save(sprintf('/procdata/parksh/_marmoset/behavior/eyeData_Ca/pupilSize_BPM_%s.mat', nameSubj), 'Results_session')

end

% histogram(catPupilMean(trial_id_marmosetFace))

%% possible figures
% pupil size to 5 different categories of images for two animals
nameSubj = 'Tabla'; %'Max'; %'Tabla';
load(sprintf('/procdata/parksh/_marmoset/behavior/eyeData_Ca/pupilSize_BPM_%s.mat', nameSubj))

% Results_session(1).condName_all(Results_session(1).conditionType) % this prints out the name of conditions in this particular session

catPupilMean_session = cat(2, Results_session.matPupilMean_Cond);

for iCondition = 1:10
    pupilMean(iCondition).matTrials = cat(1, catPupilMean_session{iCondition,:});
    pupilMean(iCondition).mean = mean(pupilMean(iCondition).matTrials);
    pupilMean(iCondition).ste = std(pupilMean(iCondition).matTrials)./sqrt(length(pupilMean(iCondition).matTrials)-1);
end
% Need to check the condition #3, 7, 8 for differences across daily
% sessions that may require further normalization and be cautious to
% interpret the average & STE for those conditions

catMean = cat(1, pupilMean([1 2 5 6]).mean);
catSTE = cat(1, pupilMean([1 2 5 6]).ste);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
hold on
line([1:4; 1:4], [catMean+catSTE, catMean-catSTE]', 'Color', 'b', 'LineWidth', 2)
plot(catMean, 'bo-', 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', 'w')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', 1:4,  'XTickLabel', {'HF', 'MF', 'NF-Obj', 'F-Obj'})
title(sprintf('%s: Average Pupil Size', nameSubj))















