% demo_pupilSizeBPM_20210422_SHP.m
%
% To compute pupil size during stimulus presentation for each trial
%

%% Define filename.
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
setSubjName = {'Tabla', 'Max'};
for iSubj = 1:length(setSubjName)
    clear Results*
    
    nameSubj =  setSubjName{iSubj}; %'Max'; %'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
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
        % checking valid trials: first, based on "TrialError" parameter
        locValidTrials = find(cat(1, data.TrialError)>1); 
        if sum(ismember(locValidTrials, 1))>0 % excluding the first trial becuase it's repeated as the last trial
            locValidTrials = locValidTrials(~ismember(locValidTrials,1));
        end
        
        clear results*
        results_time = NaN(length(locValidTrials), 2);

        setValidTrials=[];
        iValidTrial = 0;
        for iTT = 1:length(locValidTrials) % which trial?
            % while 1
            
            %     iTT = iTT+1;
            iTrial = locValidTrials(iTT);
            
            % retrieve pupil data of this trial
            pData = data(iTrial).AnalogData.General.Gen1; % pupil data
            
            % get the time information from BehavioralCodes
            % Note for CodeNumbers: 20 (stimulus on), 55 (stimulus off)
            index_stimOn = find(data(iTrial).BehavioralCodes.CodeNumbers == 20);
            index_stimOff = find(data(iTrial).BehavioralCodes.CodeNumbers == 55);
            
            time_stimOn = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOn));
            time_stimOff = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOff));
            
            % Only include the trials with decreasing pupil size change
            % Exclude trials with unusual pupil size change using a simple
            % polynomial fit for now
            tempEye = pData(time_stimOn:time_stimOff);
%             tempEye(tempEye<-4.5) = NaN; % invalid points
            baseline = mean(pData(time_stimOn:time_stimOn+150));            
            x = [1:size(tempEye,1)]';
            p = polyfit(x, tempEye-baseline, 1);
            f = polyval(p, x);
            
            figure(101); set(gcf, 'Color', 'w'); clf;
                plot(tempEye-baseline);
                hold on;
                plot(x, f, 'r');
                legend('pupil', 'polyfit')
                title('Pupil size change during stimulus presentation')
                xlabel('Time (ms)')
                
                input('')
            
            if p(1) < -0.0001 && mean(tempEye(400:500)-baseline)<0 % include trials with decreasing trend of pupil size change (constriction upon stimulus presentation)                        
                
                tLoc_loss = tempEye<-4.5;
                if sum(tLoc_loss)>0
                    continue;
                else
                
                iValidTrial = iValidTrial+1;
                    
                setValidTrials = cat(1, setValidTrials, iTrial);
                % stack the time information for each trial
                results_time(iTT, 1) = time_stimOn;
                results_time(iTT, 2) = time_stimOff;
                
                % get pupil size during stimulus presentation (between stim on and stim off timing)
                results_pupil(iValidTrial).pupilData_raw = tempEye; %pData(time_stimOn:time_stimOff);
                results_pupil(iValidTrial).pupilData_baseline = baseline; %mean(pData(time_stimOn:time_stimOn+150));
                results_pupil(iValidTrial).pupilData_baselineSubtracted = tempEye-baseline; %results_pupil(iTT).pupilData_raw - results_pupil(iTT).pupilData_baseline; % pData(time_stimOn:time_stimOff);
                
%                 sortedPupil = sort(results_pupil(iTT).pupilData_norm);
%                 results_pupil(iTT).pupilData_min150ms_mean = mean(sortedPupil(1:150));
                
%                 figure(101); set(gcf, 'Color', 'w'); clf;
%                 plot(results_pupil(iTT).pupilData_baselineSubtracted);
%                 hold on;
%                 plot(x, f, 'r');
%                 legend('pupil', 'polyfit')
%                 title('Pupil size change during stimulus presentation')
%                 xlabel('Time (ms)')
% 
%                 input('')
                end

            else % exclude any trial with increasing trend of pupil size change
                results_time(iTT, 1) = NaN;
                results_time(iTT, 2) = NaN;
                continue;
            end            
            
        end
        
        if ~exist('results_pupil', 'var')
%             Results_trial.setValidTrials = setValidTrials;
            continue;
        end
        
        %% Normalize the pupil change further using the maximum constriction across all trials
        temp = cat(1, results_pupil.pupilData_baselineSubtracted);
        sessionMaxConstriction = min(temp);

        for ii = 1:length(results_pupil)
            results_pupil(ii).pupilData_norm_relativeConstriction= (results_pupil(ii).pupilData_baselineSubtracted)./sessionMaxConstriction;
            sortedPupil = sort(results_pupil(ii).pupilData_norm_relativeConstriction, 'descend');
            results_pupil(ii).pupilData_norm_relConstriction_median = median(sortedPupil(1:150));
        end
        
        % concatenate the pupil size summary statistics from each trial across trials
        catPupilMean=[]; catPupilSte=[]; 
        catPupilMean = cat(1, results_pupil.pupilData_norm_relConstriction_median); % average (mean) pupil size
%         catPupilSte = cat(1, results_pupil.pupilData_ste); % standard error of the mean

        
        % % quick histogram
        % figure;
        % histogram(catPupilMean)
        % histogram(catPupilMean, 30)
        
        Results_trial.setValidTrials = setValidTrials;
        Results_trial.catPupilMean = catPupilMean;
%         Results_trial.catPupilSte = catPupilSte;
        Results_trial.results_pupil = results_pupil;
        Results_trial.results_time = results_time;
        Results_trial.info = ['results_pupil contains raw pupil size data, average and standard error for each trial.', 'results_time contains stimulus on (1st column) and off (2nd column) time stamps for each trial'];
        
        %% Condition information
        trial_conditions=[];
        trial_conditions = cat(1, data(setValidTrials).Condition); %cat(1, data(locValidTrials).Condition);
        
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
        
        save(sprintf('/procdata/parksh/_marmoset/behavior/eyeData_Ca/pupilSizeNorm_BPM_%s.mat', nameSubj), 'Results_session')
        
        fprintf(1, '%s, file #%d/%d: # of valid trials: %d: ...saved \n', nameSubj, iFile, length(setFilename), length(setValidTrials))
%         end
    end
    
end
% histogram(catPupilMean(trial_id_marmosetFace))

%% possible figures
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs';
% pupil size to 5 different categories of images for two animals
nameSubj = 'Tabla'; %'Max'; %'Tabla';
load(sprintf('/procdata/parksh/_marmoset/behavior/eyeData_Ca/pupilSizeNorm_BPM_%s.mat', nameSubj))

% Results_session(1).condName_all(Results_session(1).conditionType) % this prints out the name of conditions in this particular session

catPupilMean_session = cat(2, Results_session.matPupilMean_Cond);

for iCondition = 1:10
    pupilMean(iCondition).matTrials = cat(1, catPupilMean_session{iCondition,:});
    pupilMean(iCondition).mean = mean(pupilMean(iCondition).matTrials);
    pupilMean(iCondition).median = median(pupilMean(iCondition).matTrials);
    pupilMean(iCondition).ste = std(pupilMean(iCondition).matTrials)./sqrt(length(pupilMean(iCondition).matTrials)-1);
end
% Need to check the condition #3, 7, 8 for differences across daily
% sessions that may require further normalization and be cautious to
% interpret the average & STE for those conditions

catMean = cat(1, pupilMean([1 2 5 6]).median);
catSTE = cat(1, pupilMean([1 2 5 6]).ste);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
hold on
line([1:4; 1:4], [catMean+catSTE, catMean-catSTE]', 'Color', 'b', 'LineWidth', 2)
plot(catMean, 'bo-', 'LineWidth', 3, 'MarkerSize', 14, 'MarkerFaceColor', 'w')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', 1:4,  'XTickLabel', {'HF', 'MF', 'NF-Obj', 'F-Obj'})
xlim([0.5 4.5])
set(gca, 'LineWidth', 2, 'FontSize', 15)

ylim([0.35 0.65])%Tabla
set(gca, 'YTick', 0.35:0.1:0.65)%Tabla
% ylim([0.25 0.45])%Max
% set(gca, 'YTick', 0.25:0.1:0.45)%Max
print(gcf, fullfile(dirFig, sprintf('RelConstriction_BPM_HF_MF_NFO_FO_median_%s', nameSubj)), '-depsc')
% title(sprintf('%s: Average Pupil Size', nameSubj))

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
hold on
for iCond = 1:10
if isempty(pupilMean(iCond).matTrials)
continue;
end
plot(iCond, pupilMean(iCond).matTrials, 'bo'); 
% line([iCond+0.2 iCond+0.2], [pupilMean(iCond).mean+pupilMean(iCond).ste, pupilMean(iCond).mean-pupilMean(iCond).ste], 'Color', 'b', 'LineWidth', 2)
plot(iCond+0.2, pupilMean(iCond).median, 'ko', 'MarkerFaceColor', 'w', 'LineWidth', 2, 'MarkerSize', 12)
hold on
end
xlim([0 11])


%% Example trial figure
filename = '191113_Tabla_Ca_BPM_123049.bhv2';
dateSession = filename(1:6);
if str2num(dateSession) < 191121
dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
else
dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
end
data = mlread(fullfile(dirBHV, filename)); % mlread(filename);
iTrial = 8;
pData = data(iTrial).AnalogData.General.Gen1; % pupil data
index_stimOn = find(data(iTrial).BehavioralCodes.CodeNumbers == 20);
index_stimOff = find(data(iTrial).BehavioralCodes.CodeNumbers == 55);
time_stimOn = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOn));
time_stimOff = floor(data(iTrial).BehavioralCodes.CodeTimes(index_stimOff));
xAxis = [1:length(pData)]-time_stimOn;

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(xAxis, pData, 'k-', 'LineWidth', 2)
axis tight
ylim([-5.5 -2])
hold on
line([0 time_stimOff-time_stimOn; 0 time_stimOff-time_stimOn], repmat(get(gca, 'YLim')', 1, 2), 'Color', ones(1, 3).*0.6);
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12)

%% zoomed in of the "stimulus on" time window of the same trial above
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(pData(time_stimOn:time_stimOff), 'k-', 'LineWidth', 2)
axis tight
ylim([-4.5 -2])
hold on
line([150 150], get(gca, 'YLim'), 'Color', 'b');
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12)












