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

% flagSaveFile = 0; %1;

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %2; %1; %2; %1; %2; %1; %2; %1;

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
load(fname_caTSFOV, 'cellTS')
cellTS_BPM = cellTS;
clear cellTS

% DFL data to select the cells
fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_DFL, 'cellTS', 'cellPix')
cellTS_DFL = cellTS;
clear cellTS
indCellValid = find(cat(1, cellTS_DFL.nTrial1_total)>8); %

%%
% for iS = 1:length(setDateSession)
%     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
%     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%     load(fullfile(dirProcdata_session, 'BPM_ts_tML'));
%     
%     resultsBPM(iS).stimTiming_BPM = stimTiming_BPM;
% %     resultsBPM(iS).tS_session = tS_session;
% %     resultsBPM(iS).tS_session_stim = tS_session_stim;
% 
% %     tS_session.tS_trial(iCell, iTrial)
% %     tS_sesion_stim(iCell, iStim) %% cells x 25 conditions
% end

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

% iC = 1; %:length(indCellValid)
% iCell = indCellValid(iC);
% 
% tempCat = cat(1, cellTS_BPM(iCell,:).matTS_norm_avg);
% matTS_norm_cat = reshape(tempCat', [46 5 5]);
% matTS_norm_cat_avg = squeeze(mean(matTS_norm_cat, 2));
% 
% cMap_img = cool(25);
% figure;
% set(gcf, 'Color', 'w', 'Position', [200 1000 1270 370])
% for iCat = 1:5
%     sp(iCat) = subplot(1,6,iCat);
%     line([10 20; 10 20], [-0.5 -0.5; 2.5 2.5], 'Color', [1 1 1].*0.5);
%     hold on;
%     plot(squeeze(matTS_norm_cat(:,:,iCat)))
%     set(gca, 'colororder', cMap_img((iCat-1)*5+1:iCat*5, :))
% end
% 
% sp(6) = subplot(1,6,6); % average
% line([10 20; 10 20], [-0.5 -0.5; 2.5 2.5], 'Color', [1 1 1].*0.5);
% hold on;
% plot(matTS_norm_cat_avg)
% set(sp(6), 'colororder', cool(5))
% L = legend('HF', 'MF', 'FO', 'NFO', 'RDM');
% set(L, 'Position', [0.9133 0.7291 0.0586 0.1956]); %, 'Location', 'southeast');
% set(sp(1:5), 'YLim', [-0.5 2.5])
% set(sp(:), 'XTick', 10:10:40, 'XTickLabel', 0:1:4, 'TickDir', 'out', 'Box', 'off')


%% Check BPM results in combination with movie responses
% fig_snapshot = figure;
% set(fig_snapshot, 'Color', 'w', 'Position', [200 1000 1120 700])
for iC = 1:length(indCellValid)
    
    fig_snapshot = figure;
    set(fig_snapshot, 'Color', 'w', 'Position', [200 1000 1120 700])
    
    iCell = indCellValid(iC);
    
    tempCat = cat(1, cellTS_BPM(iCell,:).matTS_norm_avg);
    matTS_norm_cat = reshape(tempCat', [46 5 5]);
    matTS_norm_cat_avg = squeeze(mean(matTS_norm_cat, 2));
    
    % FOV
    sp(1) = subplot(2,6,[1 2]); cla;
    plot(cellPix(iCell).contourCell{1}(1,1:end), cellPix(iCell).contourCell{1}(2,1:end), ...
        'Color', 'k', 'linewidth', 1); hold on;
    text(cellPix(iCell).contourCell{1}(1,end)+5, cellPix(iCell).contourCell{1}(2,end), num2str(iCell), ...
        'color', 'k')
    set(sp(1), 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])
    title(sp(1), sprintf('%s FOV%d', nameSubj, FOV_ID))
    
    % movie 1 response
    sp(2) = subplot(2,6,3:6);
    imagesc(zscore(cellTS_DFL(iCell).matTS_movie1, 0, 2));
    set(sp(2), 'TickDir', 'out', 'XTickLabel', 20:20:120, 'box', 'off')
    title(sprintf('Cell %d/%d: Cell ID %d Movie 1', iC, length(indCellValid), iCell))
    colormap(hot)
    xlabel('Time (s)')
    
    % image response: averaged across trials for each image
    cMap_img = cool(25);
    for iCat = 1:5
        sp(6+iCat) = subplot(2,6,6+iCat); cla;
        line([10 20; 10 20], [-0.5 -0.5; 2.5 2.5], 'Color', [1 1 1].*0.5);
        hold on;
        plot(squeeze(matTS_norm_cat(:,:,iCat)));
        [a,i] = max(squeeze(matTS_norm_cat(:,:,iCat)));
        text(i, a, cellstr(string(cat(1,cellTS_BPM(iCell,(iCat-1)*5+1:iCat*5).idStim))));
        set(gca, 'colororder', cMap_img((iCat-1)*5+1:iCat*5, :))
    end
    
    % image response: averaged across trials for each category
    sp(12) = subplot(2,6,12); cla; % average of BPM
    plot(matTS_norm_cat_avg)
    hold on;
    line([10 20; 10 20], repmat(get(gca, 'YLim')', 1, 2), 'Color', [1 1 1].*0.5);
    set(sp(12), 'colororder', cool(5))
    L = legend('HF', 'MF', 'NFO', 'FO', 'RDM');
    set(L, 'Position', [0.8934 0.3548 0.0670 0.1014]); %, 'Location', 'southeast');
    
    set(sp(7:11), 'YLim', [-0.5 2.5])
    set(sp(7:12), 'XTick', 10:10:40, 'XTickLabel', 0:1:4, 'TickDir', 'out', 'Box', 'off')
    xlabel(sp(9), 'Time (s)')
    
%     input('')
%     set(sp(7:12), 'NextPlot', 'replace')   
end

flagSavePPTX = 0; %1;
if flagSavePPTX
    % save figures
    addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
    addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
    
    fname_pptx = sprintf('%s_FOV%d_acrossSessionFOV_BPMDFLresponse', nameSubj, FOV_ID); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
    exportFigsToPPTX_SHP(fname_pptx);
    
    switch lower(nameSubj)
        case 'tabla'
            dest = '/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
        case 'max'
            dest = '/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
    end
    movefile(sprintf('./%s*.pptx', nameSubj), dest);

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



%%
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
