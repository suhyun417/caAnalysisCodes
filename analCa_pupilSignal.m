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
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end


setNameSubj = {'Tabla', 'Max'};

dirSave{1} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Tabla/FOV1');
dirSave{2} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Max/FOV3');

nameSubj =  'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
% dateSession = '20191113'; %'20191125';  
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% for iSession = 1:nSession
    iSession = 12;
    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    if str2num(dateSession) < 20191121
        dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
    else
        dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
    end
    dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
    
    
    % load and plot the cell time series to check
    load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), ...
        'tS_session_stim', 'stimTiming_BPM');
    iCell = 82;
    iStim = 23;
    
    cMap_trial = hsv(size(tS_session_stim(iCell, iStim).matTS_norm, 2));
    fig_ts = figure;
    set(gca, 'ColorOrder', cMap_trial);
    hold on;
    plot([-1:0.1:3.5], tS_session_stim(iCell, iStim).matTS_norm)
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
        data(indTrial_bhv).BehavioralCodes.CodeNumbers == 20
    
    
    end
    fig_pupil = figure;
    set(gca, 'ColorOrder', cMap_trial)
    hold on;
    for iT = 1:length(tempPupilSig)
        figure(fig_pupil);
        plot(tempPupilSig{iT})
        hold on;
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
        



