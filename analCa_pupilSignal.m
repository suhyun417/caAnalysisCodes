% analCa_pupilSignal.m
%
% quick and dirty stuff to look at the pupil signal and absence of
% responses 

clear all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
    dirRawdata_archive = '/Volumes/archive_rawdata1/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
    dirRawdata_archive = '/rawdata_archive1/parksh';
end


setNameSubj = {'Tabla', 'Max'};

dirSave{1} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Tabla/FOV1');
dirSave{2} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Max/FOV3');

nameSubj =  'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
% dateSession = '20191113'; %'20191125';  
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% for iSession = 1:nSession
    iSession = 1; %3;
    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    if str2num(dateSession) < 20191121
        dirBHV = [dirRawdata_archive, '/behavior/MonkeyLogic_Ca/']; %
    else
        dirBHV = [dirRawdata, '/behavior/MonkeyLogic_Ca/']; %
    end
    dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
    
    % get the info for concatenated runs & exp types
    load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc/ConcatRuns_BPM_DFL.mat', ...
        nameSubj, dateSession)), 'paramConcat');
    infoSession = paramConcat.infoSession;
    catExpName = cellstr(cat(1, infoSession.ExpName));
    indBPM = find(contains(catExpName, 'BPM')>0);
       
    % load and plot the cell time series to check
    load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), ...
        'stimTiming_BPM', 'tS_session')
%     indTrial_org = cat(2, tS_session(1).idRunTrial(:,1), cat(1, stimTiming_BPM.indValidTrial_orgBhv));
    tS_trial_session = tS_session(1).tS_trial;
    
    for iTrial = 1:size(tS_trial_session, 2)
    
        % population response
        tempTS = cat(2, tS_trial_session(:,iTrial).matTS_norm);
        
        indRun = tS_session(1).idRunTrial(iTrial,1); %indTrial_org(iTrial,1);
%         indTrial = indTrial_org(iTrial,2);
        clear data tempEye
        
        % retrieve and quantify eye signal
        tempEye = stimTiming_BPM(indRun).analog.eye(round(stimTiming_BPM(indRun).t_org.stimOnset(iTrial)):round(stimTiming_BPM(indRun).t_org.blankOnset_afterStim(iTrial)), :);
        tempEye(tempEye<-5) = NaN;
        
        propLoss = sum(isnan(tempEye))./size(tempEye,1);
        x = [1:size(tempEye,1)]';
        p = polyfit(x, tempEye(:,3), 1);
        f = polyval(p, x);        

        %         fname_bhv = fullfile(dirBHV, [infoSession(indBPM(indRun)).MLFilename '.bhv2']);
%         [data] = mlread(fname_bhv);
%         stimOn = round(data(indTrial).BehavioralCodes.CodeTimes(data(indTrial).BehavioralCodes.CodeNumbers == 20));
%         stimOff = round(data(indTrial).BehavioralCodes.CodeTimes(data(indTrial).BehavioralCodes.CodeNumbers == 55));
%         tempEye = data(indTrial).AnalogData.Eye(stimOn:stimOff, :);
%         tempPS = data(indTrial).AnalogData.General.Gen1(stimOn:stimOff); %
%         tempPS(tempPS<-5) = NaN;
        
        figure(101); set(gcf, 'Color', 'w'); clf;
        sp(1) = subplot('Position', [0.05 0.6 0.5 0.35]);
        plot(tempEye(:, 1:2))
        legend('x', 'y')
        title(sp(1), 'Eye gaze during stimulus presentation')
        
        sp(2) = subplot('Position', [0.05 0.15 0.5 0.35]);
        plot(tempEye(:, 3));
        hold on;
        plot(x, f, 'r');
        legend('pupil', 'polyfit')
        title(sp(2), 'Pupil size change during stimulus presentation')
        xlabel('Time (ms)')
        
        sp(3) = subplot('Position', [0.65 0.15 0.3 0.8]);
        plot([-1:0.1:3.5], tempTS);
        xlabel('Time from stimulus onset (s)')
        ylabel('Norm. resp. (z)')
        title(sp(3), sprintf('%s %s: Trial #%d: Stimulus ID %d', nameSubj, dateSession, iTrial, tS_trial_session(1,iTrial).idStim))
        
        axis(sp(:), 'tight');
        text(sp(2), 1, max(get(sp(2), 'YLim')), sprintf('slope = %2.3f', p(1)), 'VerticalAlignment', 'top');
        
        input('')
        
    end
    
    
    
    
    
            x = [1:length(y)]';
        p = polyfit(x, y, 1);
    
    tempPupilSig = {};
    for iT = 1:length(tS_session_stim(iCell, iStim).indTrial_org)
        indRun = tS_session_stim(iCell, iStim).indTrial_org(iT,1);
        indTrial_bhv = stimTiming_BPM(indRun).indValidTrial_orgBhv(tS_session_stim(iCell, iStim).indTrial_org(iT,2));
        
        fname_bhv = fullfile(dirBHV, [infoSession(indBPM(indRun)).MLFilename '.bhv2']);
        [data] = mlread(fname_bhv);
        tempPupilSig{iT} = data(indTrial_bhv).AnalogData.General.Gen1; %data(tS_session_stim(iCell, iStim).indTrial_org(iT,2)).AnalogData.General.Gen1;
        
        % CodeNumbers: 20 for Stim Onset, 25 for Blank After Stim Off
        stimOn = round(data(indTrial_bhv).BehavioralCodes.CodeTimes(data(indTrial_bhv).BehavioralCodes.CodeNumbers == 20));
        stimOff = round(data(indTrial_bhv).BehavioralCodes.CodeTimes(data(indTrial_bhv).BehavioralCodes.CodeNumbers == 55));
        tempPupilSig_stimOn{iT} = tempPupilSig{iT}(stimOn:stimOff);
    end
    
    
    
    
    %%
    % load and plot the cell time series to check
    load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), ...
        'tS_session_stim', 'stimTiming_BPM');
    iCell = 82; %21; %82;
    iStim = 23; %13; %23;
    
    setMarkers = {'o', '^', 'sq', 'x', '+'};     
    cMap_trial = hsv(size(tS_session_stim(iCell, iStim).matTS_norm, 2));
    fig_ts = figure;
%     set(groot, 'DefaultAxesColorOrder', hsv(4), 'DefaultAxesLineStyleOrder', '-|--|:');
    set(gca, 'ColorOrder', cMap_trial);
    hold on;
    for iTrial = 1:size(tS_session_stim(iCell, iStim).matTS_norm, 2)
        figure(fig_ts);
        plot([-1:0.1:3.5], tS_session_stim(iCell, iStim).matTS_norm(:, iTrial), 'o-', 'Marker', setMarkers{mod(iTrial, 5)+1}, 'Color', cMap_trial(iTrial,:));
        hold on
    end
%     plot([-1:0.1:3.5], tS_session_stim(iCell, iStim).matTS_norm)
    title(sprintf('%s %s: Cell ID #%d, Stimulus ID #%d', nameSubj, dateSession, iCell, iStim))
    
    
    % get the info for concatenated runs & exp types
    load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL.mat'), 'paramConcat');
    
    infoSession = paramConcat.infoSession;
    catExpName = cellstr(cat(1, infoSession.ExpName));
    indBPM = find(contains(catExpName, 'BPM')>0);
    
%     stimTiming_BPM(tS_session_stim(iCell, iStim).indTrial_org(iT,1)).indValidTrial_orgBhv(tS_session_stim(iCell, iStim).indTrial_org(iT,2));
    
    tempPupilSig = {};
    for iT = 1:length(tS_session_stim(iCell, iStim).indTrial_org)
        indRun = tS_session_stim(iCell, iStim).indTrial_org(iT,1);
        indTrial_bhv = stimTiming_BPM(indRun).indValidTrial_orgBhv(tS_session_stim(iCell, iStim).indTrial_org(iT,2));
        
        fname_bhv = fullfile(dirBHV, [infoSession(indBPM(indRun)).MLFilename '.bhv2']);
        [data] = mlread(fname_bhv);
        tempPupilSig{iT} = data(indTrial_bhv).AnalogData.General.Gen1; %data(tS_session_stim(iCell, iStim).indTrial_org(iT,2)).AnalogData.General.Gen1;
    
        % CodeNumbers: 20 for Stim Onset, 25 for Blank After Stim Off
        stimOn = round(data(indTrial_bhv).BehavioralCodes.CodeTimes(data(indTrial_bhv).BehavioralCodes.CodeNumbers == 20));
        stimOff = round(data(indTrial_bhv).BehavioralCodes.CodeTimes(data(indTrial_bhv).BehavioralCodes.CodeNumbers == 55));
        tempPupilSig_stimOn{iT} = tempPupilSig{iT}(stimOn:stimOff);
    end
    
        setMarkers = {'o', '^', 'sq', 'x', '+'};     
    cMap_trial = hsv(size(tS_session_stim(iCell, iStim).matTS_norm, 2));
    fig_pupil = figure;
    set(gca, 'ColorOrder', cMap_trial)
    hold on;
    for iT = 1:length(tempPupilSig)
        figure(fig_pupil);
        plot(tempPupilSig{iT}, 'o-', 'Marker', setMarkers{mod(iT, 5)+1}, 'Color', cMap_trial(iT,:));
        hold on;
    end
    
    for iT = 1:length(tempPupilSig_stimOn)
        y = tempPupilSig_stimOn{iT};
        x = [1:length(y)]';
        p = polyfit(x, y, 1);
        matSlope(iT, :) = p;
    end    
    setValidTrial = find(matSlope(:,1)<-0.0001);
    
    setMarkers = {'o', '^', 'sq', 'x', '+'};
    cMap_trial = hsv(size(tS_session_stim(iCell, iStim).matTS_norm, 2));
    fig_ts = figure;
%     set(groot, 'DefaultAxesColorOrder', hsv(4), 'DefaultAxesLineStyleOrder', '-|--|:');
    set(gca, 'ColorOrder', cMap_trial);
    hold on;
    for iTrial = 1:length(setValidTrial) %size(tS_session_stim(iCell, iStim).matTS_norm, 2)
        idTrial = setValidTrial(iTrial);
        figure(fig_ts);
        plot([-1:0.1:3.5], tS_session_stim(iCell, iStim).matTS_norm(:, idTrial), 'o-', 'Marker', setMarkers{mod(iTrial, 5)+1}, 'Color', cMap_trial(iTrial,:));
        hold on
    end
    
    % get the info for concatenated runs & exp types
    load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL.mat'), 'paramConcat');
    
    infoSession = paramConcat.infoSession;
    catExpName = cellstr(cat(1, infoSession.ExpName));
    indBPM = find(contains(catExpName, 'BPM')>0);
       
    for iRun = 1:length(indBPM)
        
    fname_bhv = fullfile(dirBHV, [infoSession(indBPM(iRun)).MLFilename '.bhv2']);
    [data, MLConfig, TrialRecord] = mlread(fname_bhv);    
    
    figure(200);
    set(gcf, 'Color', 'w')
    for iTrial = 1:length(data)
        plot(data(iTrial).AnalogData.Eye); hold on;
        plot(data(iTrial).AnalogData.General.Gen1);
        hold on
        line(repmat(data(iTrial).BehavioralCodes.CodeTimes, 1, 2)', repmat(get(gca, 'YLim')', 1, length(data(iTrial).BehavioralCodes.CodeTimes)))
        title(sprintf('%s %s: Run #%d/%d: Trial #%d/%d', nameSubj, dateSession, iRun, length(indBPM), iTrial, length(data)));
        input('')
        clf;
    end
    end
    
    %% DFL 
    infoSession = paramConcat.infoSession;
    catExpName = cellstr(cat(1, infoSession.ExpName));
    indDFL = find(contains(catExpName, 'DFL')>0);
       
    for iRun = 1:length(indDFL)
        
    fname_bhv = fullfile(dirBHV, [infoSession(indDFL(iRun)).MLFilename '.bhv2']);
    [data, MLConfig, TrialRecord] = mlread(fname_bhv);
    
    
    figure(100);clf;
    set(gcf, 'Color', 'w')
    for iTrial = 1:length(data)
        plot(data(iTrial).AnalogData.Eye); hold on;
        plot(data(iTrial).AnalogData.General.Gen1);
        hold on
        line(repmat(data(iTrial).BehavioralCodes.CodeTimes, 1, 2)', repmat(get(gca, 'YLim')', 1, length(data(iTrial).BehavioralCodes.CodeTimes)))
        title(sprintf('%s %s: Run #%d/%d: Trial #%d/%d', nameSubj, dateSession, iRun, length(indBPM), iTrial, length(data)));
%         input('')
%         clf;
    end
    figure;
    subplot(2,1,1)
    plot(data(iTrial).AnalogData.Eye(5000:20000,1))
    subplot(2,1,2)
    plot(data(iTrial).AnalogData.General.Gen1(5000:20000));
    end
    
                
        %% Read source data and compute center coordinates of cells
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        validIndCell = [];
        validIndCell(:,1) = 1:length(neuron.ids);
        if strcmpi(nameSubj, 'max')
            load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
            validIndCell = indCell.validCell;
        end
        
        [center] = neuron.estCenter();
        center = center(validIndCell, :);
        




