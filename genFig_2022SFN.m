% genFig_2022SFN.m
%
% 2022/09/26 SHP
%   - generate figures for 2022 SFN talk
%   - part of this script is to generate a figure for 2022 NIMH Training Day
%  poster
% 

clear all;

%% Settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
%     addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/NIFVAULT/projects';
        directory.procdata = '/Volumes/NIFVAULT/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/NIFVAULT/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/projects';
        directory.procdata = '/nifvault/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/nifvault/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end
dirFig = '/nifvault/projects/parksh/0Marmoset/Ca/_labNote/_figs';

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
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')

fname_caTSFOV_RS = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS_RS', 'resultsRS')

%%
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 
% indCellValid_trial = find(cat(1, cellTS.nTrial1_total)>8); % cells that have more than 8 trials for movie 1
% 
for ii = 1:length(cellTS)
    minsnrmovie1(ii, 1) = min(cellTS(ii).snr_movie1);
end

% indCellValid_snr = find(minsnrmovie1>0.1);   

% indCellValid = intersect(indCellValid_trial, indCellValid_snr);

clear matAvg*
for iCell = 1:length(indCellValid)   
    
    matAvgTS1(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie1)'; % now it's scaled dF
    matAvgTS2(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie2)'; %
    
    steAvgTS1(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie1)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie1, 1)-1))'; % now it's scaled dF
    steAvgTS2(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie2)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie2, 1)-1))'; 

end

[coeff, score, latent, tsquared, explained] = pca(zscore(matAvgTS1)');
[sortedScore, indCell] = sort(score(:,1), 'descend');
[sortedScore2, indCell2] = sort(score(:,2), 'descend');

explained(1:10)

%% Read source data
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


% get the contours and image field of view
% neuron_b = neuron.batches{1}.neuron;

% Generate cell location map within FOV
thr = 0.5; % the lower the smaller (more centralized) the contour
cellColor = [1 1 1];
widthContour = 1;
[d1,d2] = size(neuron.Cn);

figure;
imagesc(zeros(d1, d2)); % background
colormap(gray);
caxis([0 0.1]);
hold on;

CC = cell(size(neuron.A, 2),1);
CR = cell(size(neuron.A, 2),2);
for i = 1:size(neuron.A ,2)
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
title(sprintf('%s: %s', nameSubj, dateSession))
% %save
% print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');

% end
%
Coor = neuron.get_contours(thr); 
CC = Coor;
for i = 1:size(Aor,2)
    %         cont = medfilt1(Coor{i}')';
    cont = Coor{i};
    if size(cont,2) > 1
        plot(cont(1,1:end),cont(2,1:end),'Color',cmap(i+size(Aor,2),:), 'linewidth', ln_wd); hold on;
    end
end


%% Population responses with clustering
load('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/', 'DFL_TS_clustering.mat', 'Clustering_all', 'paramClustering')


iSubj = 2;

k = 5;
[a, b] = min(sum(Clustering_all(iSubj).resultKMeans(k-1).cell_sumD))
IDXdfl = Clustering_all(iSubj).resultKMeans(k-1).cell_indCluster(:, b);
[sortedIDXdfl, indCelldfl] = sort(IDXdfl);
clear indCell_sort
for iType = 1:k
indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
end
%
% Plotting
fig_summary_DFL = figure;
set(gcf,  'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1085 750])
clear sp

% 2. Clustering results on 2-d PC space
figure(fig_summary_DFL);
sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
cMap_sort = hsv(k);
for iType = 1:k
plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
hold on;
end
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Location', 'best')
xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
set(gca, 'TickDir', 'out')
box off
axis square
title('Clustering based on movie response on PC space')
% 3. Clustering results on imaging field of view
figure(fig_summary_DFL);
sp(3) = subplot('Position', [0.55 0.1 0.4 0.4]);
cMap_sort = hsv(k);
[d1 d2] = size(infoCells(1).imgFOV);
% imagesc(imgFOV);
% colormap(sp(3), gray);
% hold on;
for iType = 1:k
for iC = 1:size(indCell_sort{iType}, 1)
Coor = cellPix(indCell_sort{iType}(iC, 1)).contourCell{1};
plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iType, :)); hold on;
end
end
axis on
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20], 'Color', 'k')


figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(zscore(matAvgTS1(:, indCelldfl))')
colormap(hot)
set(gca, 'CLim', [-2 10])
set(gca, 'YTick', find(diff(sortedIDXdfl)>0), 'XTickLabel', 20:20:120, 'TickDir', 'out')

%% Resting state    
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;

fname_caTSFOV_RS = fullfile(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_RS, 'cellTS_RS', 'resultsRS')








%% 2022 Training Day figures
% %% Intro figure: three different movie-driven ts of macaque face cells
% load('/nifvault/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'matTS_FP')
% 
% setCellIDs = {'06Dav', '25Dav', '27Dav'}; % three example neurons from the face patch ML
% tLoc = find(contains(matTS_FP.catChanID, setCellIDs));  
% 
% cellTS = matTS_FP.matFR_SU_10hz(1:3000, tLoc);
% cellTS_rs = resample(cellTS, 1, 10);
% 
% % Plot
% cMap_cell = [96 216 54; 255 67 161; 0 162 253]./255; % three example neurons
% 
% fig_intro = figure;
% set(fig_intro, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1000 500 650 480]);
% % taxis = 2.4:2.4:900;
% for iCell = 1:3
%     figure(fig_intro);
%     plot(cellTS_rs(:, iCell)+1*(iCell-1), 'LineWidth', 2, 'Color', cMap_cell(iCell,:));
%     hold on;
% end
% set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
% print(fig_intro, fullfile(dirFig, sprintf('macaqueCellMovie1_%s_%s_%s', setCellIDs{1}, setCellIDs{2}, setCellIDs{3})), ...
%     '-depsc')
% 
% 
% %% Correlation values
% matR = corr(cellTS, fmriTS, 'rows', 'complete', 'type', 'spearman');
% % matR = matR*(-1); %becuase of the MION
% 
% %% face regressor ts
% taxis = 2.4:2.4:900;
% 
% faceRegressor_mion([1:7,126:126+7, 251:251+7]) = NaN;
% faceRegressor_mion_norm = (faceRegressor_mion-nanmean(faceRegressor_mion))./nanstd(faceRegressor_mion);
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 100 650 160]);
% plot(taxis, faceRegressor_mion_norm, '-', 'LineWidth', 1, 'Color', [154 205 50]./255);
% set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
% set(gca, 'YColor', 'none')
% print(gcf, fullfile(dirFig, 'multipleFP_sFig_faceRegressor_mion'), '-depsc')

