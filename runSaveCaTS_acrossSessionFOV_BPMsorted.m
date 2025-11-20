% runSaveCaTS_acrossSessionFOV_BPMsorted.m
%
% 2025/11/19
% added the baseline-normalization matrix across sessions to be accumulated
% 2023/02/24
% save the BPM time series for the cells aligned across sessions for a given FOV


clear all;


%% Directory settings
directory = setDir_shp;
dirProjects = directory.dirProjects;
dirProcdata = directory.dirProcdata;
dirRawdata = directory.dirRawdata;
dirFig = directory.dirFig;


addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');


flagSaveFile = 0; %1; %0; %1;

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %1; %2; %1;

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

%%
resultsBPM = struct([]);
for iS = 1:length(setDateSession)
    dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    load(fullfile(dirProcdata_session, 'BPM_ts_tML'));
    
%     resultsBPM(iS).tS_session = tS_session;
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
% baselineNorm_param: '20251118_baseline from the entire run is used to normalize and termed "baselineNorm"'
% matTS_baselineNorm: [46×11 double]
% matResp_baselineNorm: [11×11 double]
% matRespAvg_baselineNorm: [11×1 double]

cellTS = struct([]);
cellPix = struct([]);

for iCell = 1:size(cellIDAcrossDay, 1)
    curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
    curCells_id = cellIDAcrossDay(iCell, curCells_session);
    
    for iStim = 1:25 
        tempMatTS = []; tempMatTS_norm = []; tempMatTS_basenorm = []; indTrial_org = {}; nTrial = [];
        for iSetCell = 1:length(curCells_id)
            
            if iSubj == 1 && ismember(curCells_session(iSetCell), [11 12]) % last two sessions of Tabla's are not merged b/c of stimulus set
                continue;
            end
            
            if iSubj == 2 && ismember(curCells_session(iSetCell), [1 2]) % first two sessions of Max's are not merged b/c of 500ms stim on time
                continue;
            end
            
            
            tempMatTS = cat(1, tempMatTS, resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS');
            tempMatTS_norm = cat(1, tempMatTS_norm, resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS_norm');

            tempMatTS_basenorm = cat(1, tempMatTS_basenorm, resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS_baselineNorm');
            
            indTrial_org{iSetCell} = resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).indTrial_org;
            nTrial(iSetCell, 1) = size(indTrial_org{iSetCell}, 1);
            
            curIdStim = resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).idStim;
            
        end
        cellTS(iCell, iStim).CellIDAcrossSession = cat(2, curCells_session', curCells_id');
        cellTS(iCell, iStim).indTrial_org = indTrial_org;
        cellTS(iCell, iStim).idStim = curIdStim; %resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).idStim;
        cellTS(iCell, iStim).nTrial = nTrial;
        cellTS(iCell, iStim).matTS = tempMatTS;
        cellTS(iCell, iStim).matTS_norm = tempMatTS_norm;
        cellTS(iCell, iStim).matTS_norm_avg = mean(tempMatTS_norm);
        cellTS(iCell, iStim).matTS_norm_ste = std(tempMatTS_norm)./sqrt(size(tempMatTS_norm, 1)-1);
        cellTS(iCell, iStim).matTS_baselinenorm = tempMatTS_basenorm;
        cellTS(iCell, iStim).matTS_baselinenorm_avg = mean(tempMatTS_basenorm);
        cellTS(iCell, iStim).matTS_baselinenorm_ste = std(tempMatTS_basenorm)./sqrt(size(tempMatTS_basenorm, 1)-1);
        
    end
end

%%
if flagSaveFile
    fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    save(fname_caTSFOV, 'cellTS') %, 'cellPix')
%     fprintf(1, '\n Saving files for %s FOV %d: in total %d/%d cells aligned \n', nameSubj, FOV_ID, countValidCell, length(flagDone))
end

%temporary save (due to 
fname_caTSFOV = sprintf('/home/parks23/Research/0Marmoset/Ca/tempData/%s_FOV%d_BPMsorted_baselineNormAdded.mat', nameSubj, FOV_ID);
save(fname_caTSFOV, 'cellTS') %, 'cellPix')




% for iC = 1:length(indCellValid)
%     iCell = indCellValid(iC);
%     
%     curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
%     curCells_id = cellIDAcrossDay(iCell, curCells_session);
%     
%     % quick and dirty check of cells across session to the same stimulus
%     cMap_s = turbo(length(curCells_id));
%     figure(100);
%     set(gcf, 'Position', [150 270 1450 1000]); clf;
%     for iStim = 1:25 %size(resultsBPM(curCells_session(1)).tS_session_stim, 2)
%         sp(iStim) = subplot(5,5,iStim);
%         for iSetCell = 1:length(curCells_id)
%             plot(resultsBPM(curCells_session(iSetCell)).tS_session_stim(curCells_id(iSetCell),iStim).matTS_norm, 'Color', cMap_s(iSetCell, :))
%             hold on
%         end        
%     end
%     title(sp(1), sprintf('Cell ID:%d', iCell))
%     input('')
% end

% 
% for iCell = 1:size(cellIDAcrossDay, 1)
%     curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
%     curCells_id = cellIDAcrossDay(iCell, curCells_session);
%     
%     tempMatTS1 = []; tempMatTS2 = []; nTrial1 = []; nTrial2 = []; stackCellPix = []; contourCell = {};
%     for iSetCell = 1:length(curCells_id)        
%         tempMatTS1 = cat(1, tempMatTS1, squeeze(resultsDFL(curCells_session(iSetCell)).tS_session(1).matTS_norm(:, curCells_id(iSetCell), :))');
%         nTrial1 = cat(1, nTrial1, size(tempMatTS1, 1));
% 
%         tempMatTS2 = cat(1, tempMatTS2, squeeze(resultsDFL(curCells_session(iSetCell)).tS_session(2).matTS_norm(:, curCells_id(iSetCell), :))');
%         nTrial2 = cat(1, nTrial2, size(tempMatTS2, 1));
%         
%         % spatial component 
%         stackCellPix = cat(2, stackCellPix, cellAcrossDay(curCells_session(iSetCell)).stackCell(:,curCells_id(iSetCell)));
%         
%         curContour = infoCells(curCells_session(iSetCell)).coor_0p2{curCells_id(iSetCell)};...
%         curContour_shifted = curContour + [shifts(curCells_session(iSetCell), 2); shifts(curCells_session(iSetCell), 1)]; 
%         contourCell = cat(1, contourCell, curContour_shifted);
%     end
%     
%     cellTS(iCell).idAcrossSession = cat(2, curCells_session', curCells_id');
%     cellTS(iCell).matTS_movie1 = tempMatTS1;
%     cellTS(iCell).matTS_movie2 = tempMatTS2;
%     cellTS(iCell).nTrial1_set = nTrial1;
%     cellTS(iCell).nTrial1_total = nTrial1(end);
%     cellTS(iCell).nTrial2_set = nTrial2;
%     cellTS(iCell).nTrial2_total = nTrial2(end);
%     
%     cellPix(iCell).idAcrossSession = cat(2, curCells_session', curCells_id');
%     cellPix(iCell).stackCellPix = stackCellPix;
%     cellPix(iCell).repPix = stackCellPix(:,1);
%     cellPix(iCell).contourCell = contourCell;
% end



