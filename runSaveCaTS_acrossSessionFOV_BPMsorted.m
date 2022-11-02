% runSaveCaTS_acrossSessionFOV_BPMsorted.m
%
% 2022/11/1
% save the BPM time series for the cells aligned across sessions for a given FOV


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

iSubj = 1; %2; %1;

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
for iS = 1:length(setDateSession)
    dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    load(fullfile(dirProcdata_session, 'BPM_ts_tML'));
    
%     tS_sesion_stim(iCell, iStim) %% cells x 25 conditions
    %%%%%%%%%% bis hier %%%%%%%%%%%%%%5
    
    
%     resultsDFL(iS).tS_session = tS_session;
end

% %%
% cellTS = struct([]);
% cellPix = struct([]);
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


% %%
% if flagSaveFile
%     fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
%     save(fname_caTSFOV, 'cellTS', 'cellPix')
%     fprintf(1, '\n Saving files for %s FOV %d: in total %d/%d cells aligned \n', nameSubj, FOV_ID, countValidCell, length(flagDone))
% end
