% genFig_dataMining_BPM.m
%
% 2023/02/28 SHP
%   - playground for BPM results (flashing images)


clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/PROJECTS/parksh';
        dirProcdata = '/Volumes/NIFVAULT/PROCDATA/parksh';
        dirRawdata = '/Volumes/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

flagSaveFile = 0; %1;

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);


%% load saved files
% cell-center info pooled across days
fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_stack) %, 'cellIDAcrossDay'); %, 'stackCellCenter')

fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_cellQC, 'infoCells')

% translational shift across days
fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));  
load(fname_shifts, 'shifts')

% BPM data
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
save(fname_caTSFOV, 'cellTS')
cellTS_BPM = cellTS;
clear cellTS

% DFL data to select the cells
fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_DFL, 'cellTS')
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); %

%%
for iS = 1:length(setDateSession)
    dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    load(fullfile(dirProcdata_session, 'BPM_ts_tML'));
    
    resultsBPM(iS).stimTiming_BPM = stimTiming_BPM;
    resultsBPM(iS).tS_session = tS_session;
    resultsBPM(iS).tS_session_stim = tS_session_stim;

%     tS_session.tS_trial(iCell, iTrial)
%     tS_sesion_stim(iCell, iStim) %% cells x 25 conditions
end

% resultsBPM(1).tS_session_stim(1,1)
% idStim: 11
% indTrial: [11×1 double]
% indTrial_org: [11×2 double]
% matTS: [46×11 double]
% matTS_norm: [46×11 double]
% matAmp: [11×11 double]
% matAmp_norm: [11×11 double]
% matAvgAmp: [11×1 double]
% matAvgAmp_norm: [11×1 double]
% avgAmp: -0.0026
% avgAmp_norm: -0.2002
% matAmp_b: [10×11 double]
% matAmp_b_norm: [10×11 double]
% matAvgAmp_b: [11×1 double]
% matAvgAmp_b_norm: [11×1 double]
% avgAmp_b: -0.0864
% avgAmp_b_norm: -0.3803


for iC = 1:length(indCellValid)
    iCell = indCellValid(iC);
    
    curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
    curCells_id = cellIDAcrossDay(iCell, curCells_session);
    
    % quick and dirty check of cells across session to the same stimulus
    cMap_s = turbo(length(curCells_id));
    figure(100);
    set(gcf, 'Position', [150 270 1450 1000]); clf;
    for iStim = 1:25 %size(resultsBPM(curCells_session(1)).tS_session_stim, 2)
        sp(iStim) = subplot(5,5,iStim);
        for iSetCell = 1:length(curCells_id)
            plot(resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS_norm, 'Color', cMap_s(iSetCell, :))
            hold on
        end        
    end
    title(sp(1), sprintf('Cell ID:%d', iCell))
    input('')
end

% Check the cells across trials in one session: to check the amplitude &
% timing issue from trial to trial
iCell = 28; %90; %13;
curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
curCells_id = cellIDAcrossDay(iCell, curCells_session);
cMap_s = turbo(length(curCells_id));

iStim = 4; %3; %12; 
figure;
for iSetCell = 1:length(curCells_id)
    plot(resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS_norm, 'Color', cMap_s(iSetCell, :))
    input('')
    cla;
end

iSetCell = 5; %2; %3;
figure;
for iC = 1:size(resultsBPM(curCells_session(iSetCell)).tS_session_stim,1) % all the cells from this session
    
    [a,i] = max(resultsBPM(curCells_session(iSetCell)).tS_session_stim(iC,iStim).matTS_norm);
    plot(resultsBPM(curCells_session(iSetCell)).tS_session_stim(iC,iStim).matTS_norm);
    hold on;
    text(i, a, cellstr(string(1:length(a))));
    title(sprintf('Cell %d', iC))  
    input('')
    cla
end


%% check the eye data for certain trials
% load and plot the cell time series to check
% load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), ...
%     'stimTiming_BPM', 'tS_session')
%     indTrial_org = cat(2, tS_session(1).idRunTrial(:,1), cat(1, stimTiming_BPM.indValidTrial_orgBhv));
% tS_trial_session = tS_session(1).tS_trial;

% let's look at a particular set of trials
iCell = 13;
curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
curCells_id = cellIDAcrossDay(iCell, curCells_session);

iStim = 4; %12; %1; %12; %22; %12; 
iSession = 6; %10; %2; %4; %2; %3; %2;
setRunTrial = resultsBPM(iSession).tS_session_stim(1,iStim).indTrial_org;

for iTrial = 1:size(setRunTrial,1) %4;
idRun = setRunTrial(iTrial,1);
idTrial = setRunTrial(iTrial,2);

% retrieve and quantify eye signal
tempEye = resultsBPM(iSession).stimTiming_BPM(idRun).analog.eye(round(resultsBPM(iSession).stimTiming_BPM(idRun).t_org.stimOnset(idTrial))-1000:...
    round(resultsBPM(iSession).stimTiming_BPM(idRun).t_org.blankOnset_afterStim(idTrial))+3500, :);
tempEye(tempEye<-5) = NaN;

figure(4)
plot(tempEye); hold on;
line([1000 2000; 1000 2000], repmat(get(gca, 'YLim')', 1, 2), 'Color', 'k')
title(sprintf('Trial #%d: Run %d Trial %d', iTrial, idRun, idTrial))
input('')
end

for iTrial = 1:size(tS_trial_session, 2)
    
    % population response
    tempTS = cat(2, tS_trial_session(:,iTrial).matTS_norm);
    
    indRun = tS_session(1).idRunTrial(iTrial,1); %indTrial_org(iTrial,1);
    %         indTrial = indTrial_org(iTrial,2);
    clear data tempEye
    
    % retrieve and quantify eye signal
    tempEye = stimTiming_BPM(indRun).analog.eye(round(stimTiming_BPM(indRun).t_org.blankOnset_beforeStim(iTrial)):round(stimTiming_BPM(indRun).t_org.blankOnset_afterStim(iTrial)), :);
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

