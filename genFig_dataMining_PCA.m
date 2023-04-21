% genFig_dataMining_PCA.m
%
% 2023/03/29 SHP
% playground to compare the map of across neuron variance using PCA,
% compare the maps from different conditions
% The idea is to evaluate more continuous changes across neurons over
% cortical space than clustering or pair-wise correlations

clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/projects/parksh';
        dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
        dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
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

iSubj = 1; %2; %1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);


%% load saved files
% cells pooled across days
fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_stack, 'cellIDAcrossDay'); %, 'stackCellCenter')

% cell quality info
fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_cellQC, 'infoCells')

% translational shift across days
fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_shifts, 'shifts')

% aligned cells TS and spatial info
fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_DFL, 'cellTS', 'cellPix')
cellTS_DFL = cellTS;
clear cellTS

fname_caTSFOV_BPM = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_BPM, 'cellTS')
cellTS_BPM = cellTS;
clear cellTS

fname_caTSFOV_RS = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_RS, 'cellTS_RS') %, 'resultsRS')

%%
indCellValid = find(cat(1, cellTS_DFL.nTrial1_total)>8); %

%%
% what do I want to do?
% 1) load PCA results
% 2) map the cells based on PC score

%% Movie-based PCs
load(fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/DFL_TS_PCA.mat'))

iMovie = 1;
score = resultsPCA_DFL(iSubj).resultsPCA_run(iMovie).score;

iPC = 3;
for iPC = 1:3
% [sortedScore_abs, indCell_abs] = sort(abs(score(:,iPC)), 'descend');
[sortedScore, indCell] = sort(score(:,iPC), 'descend');


% indCell = indCell_abs;

figure;
% cMap_sort = jet(length(indCell)); %jet(length(indCell)); %hsv(k);
cMap_sort = cool(length(indCell));

imgFOV = infoCells(1).imgFOV;
[d1, d2] = size(imgFOV);
% imagesc(imgFOV); 
% hold on;
for iCell = 1:length(indCell)
    iC = indCell(iCell);
        Coor = cellPix(indCellValid(iC)).contourCell{1};
        plot(Coor(1,:)+shifts(cellPix(indCellValid(iC)).idAcrossSession(1,1),1),...
            Coor(2,:)+shifts(cellPix(indCellValid(iC)).idAcrossSession(1,1),2),...
            '.', 'Color', cMap_sort(iCell, :)); hold on;
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
title(sprintf('Ranked based on PC%d', iPC))
end


%% resting state PCs
load(fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/RS_TS_PCA.mat'), 'resultsPCA')

iSession = 1;
% for iSession = 1:12
% score = resultsPCA(iSubj).resultsPCA_session(iSession).score;
% 
for iBlock = 1:5
score = resultsPCA(iSubj).resultsPCA_2min(iSession, iBlock).score;

iPC = 1;
[sortedScore, indCell] = sort(score(:,iPC), 'descend');

figure;
% cMap_sort = jet(length(indCell)); %jet(length(indCell)); %hsv(k);
cMap_sort = cool(length(indCell));

imgFOV = infoCells(iSession).imgFOV;
[d1, d2] = size(imgFOV);
% imagesc(imgFOV); 
% hold on;
for iCell = 1:length(indCell)
    iC = indCell(iCell);
%     Coor = infoCells(iSession).coor_0p2{iC};
%     plot(Coor(1,:)+shifts(iSession,1),...
%         Coor(2,:)+shifts(iSession,2),...
%         '.', 'Color', cMap_sort(iCell, :)); hold on;
    plot(infoCells(iSession).cellCenter(iC, 2)+shifts(iSession, 1), ...
        infoCells(iSession).cellCenter(iC, 1)+shifts(iSession, 2), ...
        'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', cMap_sort(iCell, :), 'MarkerEdgeColor', 'none'); %cMap_sort(iCell, :));
    hold on;
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
title(sprintf('RS Session %d: Ranked based on PC%d', iSession, iPC))

% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_RS_PC%d_score_session%02d', nameSubj, FOV_ID, iPC, iSession)), '-depsc')

title(sprintf('RS Session %d, Block %d: Ranked based on PC%d', iSession, iBlock, iPC))
% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_RS_PC%d_score_session%02d_block%02d', nameSubj, FOV_ID, iPC, iSession, iBlock)), '-depsc')
end

%% BPM PCs
clear matAvg*
for iC = 1:length(indCellValid)
    
    clear resp* matTS_norm* tempCat
    iCell = indCellValid(iC);
    
    tempCat = cat(1, cellTS_BPM(iCell,:).matTS_norm_avg);
    matTS_norm_cat = reshape(tempCat', [46 5 5]);
    matTS_norm_cat_avg = squeeze(mean(matTS_norm_cat, 2));
    
    resp_avg_img = mean(tempCat(:, 20:30), 2)';
    resp_avg_cat = mean(matTS_norm_cat_avg(20:30, :));
    

    matAvg_img(iC, :) = resp_avg_img;
    matAvg_cat(iC, :) = resp_avg_cat;
end

[coeff_img, score_img, latent_img, tsquared_img, explained_img] = pca(matAvg_img); 
[sortedScore, indCellSorted_img] = sort(score_img(:,1), 'descend');

[coeff_cat, score_cat, latent_cat, tsquared_cat, explained_cat] = pca(matAvg_cat); 
[sortedScore, indCellSorted_cat] = sort(score_cat(:,1), 'descend');

figure;
% cMap_sort = jet(length(indCell)); %jet(length(indCell)); %hsv(k);
cMap_sort = cool(length(indCellSorted_cat));

imgFOV = infoCells(1).imgFOV;
[d1, d2] = size(imgFOV);
for iCell = 1:length(indCellSorted_cat)
    iC = indCellSorted_cat(iCell);
        Coor = cellPix(indCellValid(iC)).contourCell{1};
        plot(Coor(1,:)+shifts(cellPix(indCellValid(iC)).idAcrossSession(1,1),1),...
            Coor(2,:)+shifts(cellPix(indCellValid(iC)).idAcrossSession(1,1),2),...
            '.', 'Color', cMap_sort(iCell, :)); hold on;
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
title(sprintf('Ranked based on PC%d of categorical responses', iPC))

