% genFig_dataMining_PCA.m
%
% 2025/06/12 SHP
%  - modify the directory setting
%  - visualize PCs and PC score spatial distribution across animals
% 2023/03/29 created by SHP
% playground to compare the map of across neuron variance using PCA,
% compare the maps from different conditions
% The idea is to evaluate more continuous changes across neurons over
% cortical space than clustering or pair-wise correlations

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
% score = resultsPCA(iSubj).resultsPCA_run(iMovie).score;

for iSubj = 1:2
    score{iSubj} = resultsPCA(iSubj).resultsPCA_run(iMovie).score(:,1:3);
    coeff{iSubj} = resultsPCA(iSubj).resultsPCA_run(iMovie).coeff(:,1:3);
end

% Visualize the PC time course for both areas
figure
set(gcf, 'color', 'w', 'Position', [100 100 630 890])
SP(1) = subplot('Position', [0.1 0.7 0.85 0.25]);
SP(2) = subplot('Position', [0.1 0.4 0.85 0.25]);
SP(3) = subplot('Position', [0.1 0.1 0.85 0.25]);
plot(SP(1), [coeff{2}(:,1), coeff{1}(:,1)], 'LineWidth', 2)
plot(SP(2), [coeff{2}(:,2), coeff{1}(:,2)], 'LineWidth', 2)
plot(SP(3), [coeff{2}(:,3), coeff{1}(:,3)], 'LineWidth', 2)
cMap_areas = [102 230 179; 255 170 187]./255;
set(SP(:), 'ColorOrder', flipud(cMap_areas))
set(SP(:), 'xlim', [0 1200], 'LineWidth', 2)
set(SP(1:2), 'XTickLabel', [])
set(SP(3), 'XTickLabel', 0:20:120)
set(SP(:), 'TIckDir', 'out', 'Box', 'off')
print(fullfile(dirFig, 'DFL1_PC123_timecourse_TablaMax'), '-r300', '-dtiff')
print(fullfile(dirFig, 'DFL1_PC123_timecourse_TablaMax'), '-depsc')


% PC score spatial distribution over FOV
% make bwr colormap
n = 256;
r = [linspace(0,1,n/2), ones(1, n/2)];
g = [linspace(0,1,n/2), linspace(1,0, n/2)];
b = [ones(1,n/2), linspace(1,0,n/2)];
bwr = [r' g' b'];

% get coordinates for each point
curcenter_um=[];
for iC = 1:size(indCellValid, 1)
    curCenter = mean(cellPix(indCellValid(iC)).centerCell, 1);
    curcenter_um(iC, :) = curCenter.*3.125;
end

% Marker 'o' version
iSubj = 1;
iPC = 2;
figure;
set(gcf, 'Color', 'w')
[d1 d2] = size(infoCells(1).imgFOV);
scatter(curcenter_um(:,2), curcenter_um(:,1), 120, score{iSubj}(:,iPC), 'filled', 'o',...
    'MarkerEdgeColor', 'w', 'LineWidth', 2);
colormap(bwr)
% colorbar
axis equal
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20].*3.125, 'YLim', [0-20 d1+20].*3.125) %, 'Color', 'k', 'Box', 'off')
set(gca, 'YTick', 0:200:800)
set(gca, 'TickDir', 'out', 'Box', 'on')
if iSubj == 2
    set(gca, 'YTickLabelRotation', 270, 'XAxisLocation', 'top', 'XTickLabelRotation', 270, ...
        'XTickLabel', 600:-200:0)
end
grid on
set(gca, 'CLim', [-1 1].*23); %floor(mean(abs(round(get(gca, 'CLim'))))))
print(fullfile(dirFig, sprintf('DFL1_%s_PC%d_score_FOVum_o', nameSubj, iPC)), '-depsc')
print(fullfile(dirFig, sprintf('DFL1_%s_PC%d_score_FOVum_o', nameSubj, iPC)), '-r300', '-dtiff')

% for iC = 1:size(indCellValid, 1)
%     curCenter = mean(cellPix(indCellValid(iC)).centerCell, 1);
%     curcenter_um = curCenter*3.125;
%     plot(curcenter_um(2), curcenter_um(1), 'Marker', 'o') %, ...
% %         'MarkerFaceColor', cMap_sort(iType, :), 'MarkerEdgeColor', 'w',...
% %         'MarkerSize', 10, 'LineWidth', 2); 
%         hold on;
% end

% Marker 'X'
figure;
set(gcf, 'Color', 'w')
[d1 d2] = size(infoCells(1).imgFOV);
for iType = 1:k
    for iC = 1:size(indCell_sort{iType}, 1)
        curCenter = mean(cellPix(indCellValid(indCell_sort{iType}(iC, 1))).centerCell, 1);
        curcenter_um = curCenter*3.125;
        plot(curcenter_um(2), curcenter_um(1), 'Marker', 'x', ...
            'Color', cMap_sort(iType, :), ...
            'MarkerSize', 10, 'LineWidth', 2); hold on;
    end
end
axis equal
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20].*3.125, 'YLim', [0-20 d1+20].*3.125) %, 'Color', 'k', 'Box', 'off')
set(gca, 'YTick', 0:200:800)
set(gca, 'TickDir', 'out', 'Box', 'on')
if iSubj == 2
    set(gca, 'YTickLabelRotation', 270, 'XAxisLocation', 'top', 'XTickLabelRotation', 270, ...
        'XTickLabel', 600:-200:0)
end
grid on;
print(fullfile(dirFig, sprintf('DFL1_Clustering_%s_%dClusters_FOVum_x', nameSubj, k)), '-depsc')
print(fullfile(dirFig, sprintf('DFL1_Clustering_%s_%dClusters_FOVum_x', nameSubj, k)), '-r300', '-dtiff')



%% older version of visualization
iPC = 3;
for iPC = 1:3
% [sortedScore_abs, indCell_abs] = sort(abs(score(:,iPC)), 'descend');
[sortedScore, indCell] = sort(score(:,iPC), 'descend');
% indCell = indCell_abs;

% PC score distribution
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

