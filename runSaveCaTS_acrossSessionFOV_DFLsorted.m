% runSaveCaTS_acrossSessionFOV_DFLsorted.m
%
% 2022/10/31 SHP
% save DFL trial time series for the neurons across sessions
% using a longitudinally aligned membership of the cells


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


%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 1;

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
    load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
    
    resultsDFL(iS).tS_session = tS_session;
end

%%
cellTS = struct([]);
cellPix = struct([]);
for iCell = 1:size(cellIDAcrossDay, 1)
    curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
    curCells_id = cellIDAcrossDay(iCell, curCells_session);
    
    tempMatTS1 = []; tempMatTS2 = []; nTrial1 = []; nTrial2 = []; stackCellPix = []; contourCell = {};
    for iSetCell = 1:length(curCells_id)        
        tempMatTS1 = cat(1, tempMatTS1, squeeze(resultsDFL(curCells_session(iSetCell)).tS_session(1).matTS_norm(:, curCells_id(iSetCell), :))');
        nTrial1 = cat(1, nTrial1, size(tempMatTS1, 1));

        tempMatTS2 = cat(1, tempMatTS2, squeeze(resultsDFL(curCells_session(iSetCell)).tS_session(2).matTS_norm(:, curCells_id(iSetCell), :))');
        nTrial2 = cat(1, nTrial2, size(tempMatTS2, 1));
        
        % spatial component 
        stackCellPix = cat(2, stackCellPix, cellAcrossDay(curCells_session(iSetCell)).stackCell(:,curCells_id(iSetCell)));
        
        curContour = infoCells(curCells_session(iSetCell)).coor_0p2{curCells_id(iSetCell)};...
        curContour_shifted = curContour + [shifts(curCells_session(iSetCell), 2); shifts(curCells_session(iSetCell), 1)]; 
        contourCell = cat(1, contourCell, curContour_shifted);
    end
    
    cellTS(iCell).idAcrossSession = cat(2, curCells_session', curCells_id');
    cellTS(iCell).matTS_movie1 = tempMatTS1;
    cellTS(iCell).matTS_movie2 = tempMatTS2;
    cellTS(iCell).nTrial1_set = nTrial1;
    cellTS(iCell).nTrial1_total = nTrial1(end);
    cellTS(iCell).nTrial2_set = nTrial2;
    cellTS(iCell).nTrial2_total = nTrial2(end);
    
    cellPix(iCell).idAcrossSession = cat(2, curCells_session', curCells_id');
    cellPix(iCell).stackCellPix = stackCellPix;
    cellPix(iCell).repPix = stackCellPix(:,1);
    cellPix(iCell).contourCell = contourCell;
end

%% 
indCellValid = cat(1, cellTS.nTrial1_total)>8;

% tempATrial = cat(2, cellPix(indCellValid).repPix);
% tempATrial(~isnan(tempATrial)) = 10;
tempA = cat(2, cellPix.repPix);
tempA(~isnan(tempA)) = 1;
tempA(:, indCellValid) = tempA(:, indCellValid).*10;

imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));

figure;
imagesc(imgCells_2d)

% A_temp = full(reshape(neuron.A(:,i),d1,d2));

figure;
cmap_cell = colormap(hsv(length(cellPix)));
for iCell = 1:length(cellPix)
    plot(cellPix(iCell).contourCell{1}(1,1:end), cellPix(iCell).contourCell{1}(2,1:end), ...
        'Color', cmap_cell(iCell, :), 'linewidth', 1); hold on;
    text(cellPix(iCell).contourCell{1}(1,end), cellPix(iCell).contourCell{1}(2,end), num2str(iCell), ...
        'color', 'k')
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])

%%
if flagSaveFile
    fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    save(fname_caTSFOV, 'cellTS', 'cellPix')
    fprintf(1, '\n Saving files for %s FOV %d: in total %d/%d cells aligned \n', nameSubj, FOV_ID, countValidCell, length(flagDone))
end




%% plot
fig_supercellmovie = figure;
set(fig_supercellmovie, 'Position', [1500 1000 800 600], 'Color', 'w')
for iCell = 1:length(cellTS)
    
    figure(fig_supercellmovie);
    set(gcf, 'color', 'w')
    subplot(2,1,1)
    imagesc(cellTS(iCell).matTS_movie1);
    set(gca, 'XTickLabel', 20:20:120, 'YTick', cellTS(iCell).nTrial1_set , 'YTickLabel', setDateSession(cellTS(iCell).idAcrossSession(:,1)),'TickDir', 'out', 'Box', 'off')
    title(sprintf('Cell %d/%d: Mov 1', iCell, length(cellTS)))
    colormap(hot)
    
    subplot(2,1,2)
    imagesc(cellTS(iCell).matTS_movie2);
    set(gca, 'XTickLabel', 20:20:120, 'YTick', cellTS(iCell).nTrial2_set, 'YTickLabel', setDateSession(cellTS(iCell).idAcrossSession(:,1)),'TickDir', 'out', 'Box', 'off')
    title('Mov 2')
    xlabel('Time (s)')
    colormap(hot)
    

    input('')
end

    
    
    
    
    
    
    
    
    
