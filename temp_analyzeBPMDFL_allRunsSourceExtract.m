clear all; close all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

%% directory
nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');



%% Read source data
addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/')); 
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); 


%% Load the cell time series
load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
load(fullfile(dirProcdata_session, 'DFL_ts.mat'))


%% Stimulus timing info for BPM runs
for iRun = 1:length(tSeries_BPM)
    
    load(fullfile(dirProcdata_session, sprintf('BPM_%d_tML.mat', iRun)), 't_adj', 'infoTrial')
    
%     stimTiming_BPM(iRun).filename_org = DataFile;
    stimTiming_BPM(iRun).t_adj = t_adj;
    stimTiming_BPM(iRun).infoTrial = infoTrial;
    
%     stimTiming_BPM(iRun).indValidTrial = find(stim.trialError>1);
    
    fs = 10; %eventually should be retrieved from xml file instead of hard-coding
    locTrialStart = floor(t_adj.trialStart./(1000/fs))+1; % add one because this is frame index and the first frame is "1"
    locStimOn = floor(t_adj.stimOnset./(1000/fs))+1;
    locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./(1000/fs))+1;
    %     t_onset_adj = t_onset/1000-delay;
    %     locStimOn = floor(t_onset_adj./(1/fs));
    %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
    
    stimTiming_BPM(iRun).locCaFrame.indTrial = find(stim.trialError>1);
    stimTiming_BPM(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_BPM(iRun).locCaFrame.indTrial);
    stimTiming_BPM(iRun).locCaFrame.locStimOn = locStimOn;
    stimTiming_BPM(iRun).locCaFrame.locBlankOn_afterStim = locBlankOn_afterStim;    
    
end

%% Stimulus timing info for DFL runs
for iRun = 1:length(tSeries_DFL)
    
    load(fullfile(dirProcdata_session, sprintf('DFL_%d_tML.mat', iRun)), 't_adj', 'stim');
    
%     stimTiming_DFL(iRun).filename_org = DataFile;
    stimTiming_DFL(iRun).t_adj = t_adj;
    
    stimTiming_DFL(iRun).stim = stim;
    %         stimTiming_DFL(iRun).indValidTrial = find(stim.trialError>1);
    
    fs = 10; %eventually should be retrieved from xml file instead of hard-coding
    %         locTrialStart = floor(t_adj.trialStart./1000./(1/fs));
    locMovieOn = floor(t_adj.movieOnset./(1000/fs))+1;
    locMovieOff = floor(t_adj.sendTTL_end./(1000/fs))+1;
    locReward = floor(t_adj.reward./(1000/fs))+1;
    %         locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./1000./(1/fs));
    %     t_onset_adj = t_onset/1000-delay;
    %     locStimOn = floor(t_onset_adj./(1/fs));
    %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
    
    %         stimTiming_DFL(iRun).locCaFrame.indTrial = find(stim.trialError>1);
    %         stimTiming_DFL(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_DFL(iRun).locCaFrame.indTrial);
    stimTiming_DFL(iRun).locCaFrame.locMovieOn = locMovieOn;
    stimTiming_DFL(iRun).locCaFrame.locMovieOff = locMovieOff;
    stimTiming_DFL(iRun).locCaFrame.locReward = locReward;
        
end


%% BPM: sort the timeseries for each cell and each condition
% give IDs to the condition matrix for individual items

% temp = repmat([1:5], ncond, 1);
%     aa = cat(2, repmat(category', nDirection, 1), temp(:));

for iRun = 1:size(tSeries_BPM, 2)
    
%     indTrial = stimTiming_BPM(iRun).locCaFrame.indTrial;
    indValidTrial = find(stimTiming_BPM(iRun).locCaFrame.locTrialStart>0);
    condMat = stimTiming_BPM(iRun).stim.condMat(stimTiming_BPM(iRun).locCaFrame.indTrial(indValidTrial),:);
    
    [sortedCond, indTrialCond] = sortrows(condMat);
    setCond = unique(sortedCond(:,1));
    
    matTS = tSeries_BPM(iRun).C_raw'; %cat(1, tSeries_BPM(:,iRun).mnTS_C_raw)'; %cat(1, tSeries_BPM(:,iRun).mnTS)';
    matTS_norm = zscore(matTS);
    
    clear matTS_trial matTS_norm_trial
    for iTrial = 1:length(indValidTrial)
        matTS_trial(:,:,iTrial) = matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))-10:stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+40, :); %locEnd(iTrial), :); %locOnset_end(iTrial), :);
        matTS_norm_trial(:,:,iTrial) = matTS_norm(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))-10:stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+40, :);
    end
    
    clear resultsTrial
    for iTrial = 1:length(indValidTrial)
        matAmp = [];
        win_ms = [500 1500]; % time window of amplitude calculation (in ms: zero is stim Onset) 
        fs = 10;
        win_frame = win_ms./(1000/fs);
        matAmp = matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(1): stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(2), :)';
        matAmp_norm = matTS_norm(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(1): stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(2), :)';
                
%         resultsTrial(iTrial).cond = 
        resultsTrial(iTrial).window_ms = win_ms;
        resultsTrial(iTrial).win_frame = win_frame;
        resultsTrial(iTrial).matAmp = matAmp;
        resultsTrial(iTrial).mnAmp = mean(matAmp, 2);
        resultsTrial(iTrial).matAmp_norm = matAmp_norm;
        resultsTrial(iTrial).mnAmp_norm = mean(matAmp_norm, 2);
    end
        
    resultsBPM(iRun).resultsTrial = resultsTrial;
        %matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+50
    % matTS_trial = cat(3, sortTrial.matTS);
    % clear sortTrial
    
    for iCell = 1:size(matTS_trial, 2)
        
        for iCond = 1:length(setCond)
            tSeries_BPM_Cond(iCell, iRun, iCond).condMat = condMat(indTrialCond(sortedCond(:,1)==setCond(iCond)),:);
            tSeries_BPM_Cond(iCell, iRun, iCond).indTrial = indTrialCond(sortedCond(:,1)==setCond(iCond));
            tSeries_BPM_Cond(iCell, iRun, iCond).matTS = squeeze(matTS_trial(:, iCell, tSeries_BPM_Cond(iCell, iRun, iCond).indTrial)); % matTS_trial(:,:,sortCond(iCond).indTrial);
            tSeries_BPM_Cond(iCell, iRun, iCond).matTS_norm = squeeze(matTS_norm_trial(:, iCell, tSeries_BPM_Cond(iCell, iRun, iCond).indTrial)); % 
        end
        
    end
    
end


fig_selectivity = figure;

for iCell = 1:size(tSeries_BPM_Cond, 1)
    fprintf(1, 'Cell #%d: ', iCell)
    
%     iRun = 1;
    for iRun = 1:size(tSeries_BPM_Cond,2) 
%         curTS_catCond = [];
%         curTS_catCond = cat(2, tSeries_BPM_Cond(iCell, iRun, :).matTS);
        
        for iCond = 1:5
            figure(fig_selectivity)
            SP(5*(iRun-1)+iCond) = subplot(size(tSeries_BPM_Cond,2), 5, 5*(iRun-1)+iCond);
            plot(tSeries_BPM_Cond(iCell, iRun, iCond).matTS_norm)
%             plot(curTS_catCond(:,:,iCond))
            ylim([-3 10])
        end
        
        title(SP(3), sprintf('Cell #%d', iCell))
        
    end
    
    input('')
end


[sortedCond, indTrialCond] = sortrows(condMat);
catMatAmp = cat(2, resultsTrial(indTrialCond).matAmp);
catMnAmp = cat(2, resultsTrial(indTrialCond).mnAmp);
% figure;
% plot(catMatAmp(5,:), 'o');  
% figure;
% plot(catMatAmp(14,:), 'o');  hold on
figure(100);
for iC = 1:size(catMatAmp, 1)
    figure(100);cla;
    plot(catMatAmp(iC,:), 'o');  hold on
    line([110 209 319 429; 110 209 319 429], [-3 -3 -3 -3; 7 7 7 7], 'Color', 'm')
    title(sprintf('Cell #%d/%d:', iC, size(catMatAmp, 1)))
    input('')
end

clear tempCat
setCell = 1:size(tSeries_BPM_Cond, 1); %[5 14 37 38 46 47 48 54 76 77]; %1:size(tSeries_BPM_Cond, 1); %[5 14 37 38 46 47 48]; % [1 4 6];
for iCell = 1:length(setCell)
    idCell = setCell(iCell); %4; %1;
    for iCond = 1:5
        tempCat{iCond} = cat(2, tSeries_BPM_Cond(idCell, :, iCond).matTS);
        tempCat_norm{iCond} = cat(2, tSeries_BPM_Cond(idCell, :, iCond).matTS_norm);
    end
    for iCond = 1:5
        mnTS(:,iCond) = mean(tempCat{iCond}, 2);
        steTS(:,iCond) = std(tempCat{iCond}, [], 2)./sqrt(size(tempCat{iCond}, 2));
        
        mnTS_norm(:,iCond) = mean(tempCat_norm{iCond}, 2);
    end
    
    figure(100);
%     set(gcf, 'Position', [100 100 800 500])
    SP(1) = subplot(1,2,1);
    plot([1:51].*0.1, mnTS, 'LineWidth', 2)
    legend('HF', 'MF', 'NO', 'FO', 'RD')
    title(sprintf('Cell #%d', iCell))
    axis tight
    
    SP(2) = subplot(1,2,2);
    plot([1:51].*0.1, mnTS_norm, 'LineWidth', 2)
    legend('HF', 'MF', 'NO', 'FO', 'RD')
    title(sprintf('Cell #%d', iCell))
    ylabel('zscore')
    axis tight
    
    set(SP(:), 'Box', 'off', 'TickDIr', 'out', 'LineWidth', 2)
    
    input('')
%     print(gcf, fullfile(dirFig, sprintf('%s_%s_ConditionMerge_Cell%d', dateSession, nameSubj, iCell)), '-depsc')

end


%% DFL: sort the timeseries for each cell and each movie
catStim = cat(1, stimTiming_DFL(:).stim);
catNameMovie = {catStim.nameMovie}';
setMovie = unique(catNameMovie);

iMovie = 1;
setIndRun = find(contains(catNameMovie, setMovie{iMovie})>0);

clear matTS_movie
for iRun = 1:length(setIndRun)
    
    idRun = setIndRun(iRun);
    matTS_movie(:, :, iRun) = tSeries_DFL(idRun).mnTS_C_raw(:, stimTiming_DFL(idRun).locCaFrame.locMovieOn:stimTiming_DFL(idRun).locCaFrame.locMovieOff); 
end

indNeuron = neuron.orderROIs('snr');
setNeuron = indNeuron;
figCheck = figure;
for iCell = 1:length(setNeuron) %1:size(matTS_movie, 1)
    idCell = setNeuron(iCell);
    figure(figCheck);
    SP(1) = subplot(2, 1, 1);
    imagesc(zscore(squeeze(matTS_movie(idCell, :, :)))')
    colormap(hot)
    set(SP(1), 'CLim', [0 8]);
%     plot(squeeze(matTS_movie(iCell, :, :)))
    SP(2) = subplot(2, 1, 2);
    plot(zscore(squeeze(matTS_movie(idCell, :, :))))
    axis tight
    title(SP(1), sprintf('Cell #%d/%d (Cell ID: %d): movie %s', iCell, size(matTS_movie,1), tSeries_DFL(1).idNeuron_org(idCell), setMovie{iMovie}));

    input('')
end