% genFig_clustering.m
%
% 2022/09/26 SHP
%   - generate figures for 2022 SFN talk
%   - part of this script is to generate a figure for 2022 NIMH Training Day
%  poster
% 

clear all;

%% settings
% flagBiowulf = 0; %1; %0;
% 
% if flagBiowulf
%     directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
%     dirFig = '/data/parks20/analysis/_figs';
% else
%     ss = pwd;
%     if ~isempty(strfind(ss, 'Volume')) % if it's local
%         dirProjects = '/Volumes/VNDLab_Data/projects/parksh'; %'/Volumes/NIFVAULT/PROJECTS/parksh';
%         dirProcdata = '/Volumes/VNDLab_Data/procdata/parksh'; %'/Volumes/NIFVAULT/PROCDATA/parksh';
%         dirRawdata = '/Volumes/VNDLab_Data/rawdata/parksh'; %'/Volumes/rawdata/parksh';
%     else % on virtual machine
%         dirProjects = '/VNDLab_Data/projects/parksh'; %'/nifvault/projects/parksh';
%         dirProcdata = '/VNDLab_Data/procdata/parksh'; %'/nifvault/procdata/parksh';
%         dirRawdata = '/VNDLab_Data/rawdata/parksh'; %'/nifvault/rawdata/parksh';
%     end
% end

directory = setDir_shp;

dirProjects = directory.dirProjects;
dirProcdata = directory.dirProcdata;
dirRawdata = directory.dirRawdata;
dirFig = directory.dirFig;

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

% dirFig = '~/data/projects/parksh/0Marmoset/Ca/_labNote/_figs/'; %fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %1; %2; %1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% %% example time courses
% dateSession = setDateSession{1};
% dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% dirPreproc = fullfile(dirProcdata_session, '_preproc');
% 
% load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
% 
% pnrs = max(tSeries_DFL(1).C, [], 2)./std(tSeries_DFL(1).C_raw-tSeries_DFL(1).C, 0, 2); % peak amplitude divided by noise std
% [a, ind] = sort(pnrs, 'descend');
% 
% fig_movie = figure;
% set(fig_movie, 'Color', 'w')
% tlen = 1200;
% plot(tSeries_DFL(1).C_raw(ind(1:10), 1:tlen)'+repmat([1:10].*5, tlen, 1), 'LineWidth', 2)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
% % print(fig_movie, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_DFL1_DFLpnr10', dateSession)), '-depsc')
% 
% 
% load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
% 
% % pnrs = max(tSeries_BPM(1).C, [], 2)./std(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2); % peak amplitude divided by noise std
% % [a, ind] = sort(pnrs, 'descend');
% 
% fig_BPM = figure;
% set(fig_BPM, 'Color', 'w')
% tlen = 1200;
% plot(tSeries_BPM(1).C_raw(ind(1:10), 21:20+tlen)'+repmat([1:10].*5, tlen, 1), 'LineWidth', 2)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
% % print(fig_BPM, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_BPM1_DFLpnr10', dateSession)), '-depsc')


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
load(fname_caTSFOV_RS, 'cellTS_RS', 'resultsRS')

% clustering results
load(fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/DFL_TS_clustering.mat'), 'Clustering_all', 'paramClustering')


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

% explained(1:10)

%% spatial map
% all cells
tempA = cat(2, cellPix(:).repPix);
tempA(~isnan(tempA)) = 1;

% % selected cells
% tempA = cat(2, cellPix(indCellValid).repPix);
% tempA(~isnan(tempA)) = 1;
% % tempA(:, indCellValid) = tempA(:, indCellValid).*10;

imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));

figure;
set(gcf, 'Color', 'w')
subplot('Position', [0 0 1 1])
imagesc(imgCells_2d)
colormap(gray)
set(gca, 'CLim', [0 0.5])
truesize
box off
axis off
% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_cellMap_allCells_truesize', nameSubj, FOV_ID)), '-depsc')


%% Explained Variance plot
setK = paramClustering.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    indClust = Clustering_all(iSubj).resultKMeans(iK).cell_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
    [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(:,iK) = sum(Clustering_all(iSubj).resultKMeans(iK).cell_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering_all(iSubj).totalSS_cell; %[a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
% betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;


%% plot Clustering results
% explained variance
figure;
plot(setK, propExplained'.*100, 'ko-'); hold on
xlabel('Number of cluster (K)')
ylabel('Explained variance (%)')
set(gca, 'XTick', setK)

figure;
plot(setK(1:end-1), diff(propExplained').*100)
hold on
plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
title('difference of explained variance for each K')

% distance within each cluster
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_all(iSubj).resultKMeans(iK).cell_sumD, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')

%% Population responses with clustering: Tabla's k=5 case
% cells with more than 8 trials
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 

k = 5;
[a, b] = min(sum(Clustering_all(iSubj).resultKMeans(k-1).cell_sumD));
IDXdfl = Clustering_all(iSubj).resultKMeans(k-1).cell_indCluster(:, b);
[sortedIDXdfl, indCelldfl] = sort(IDXdfl);
clear indCell_sort
for iType = 1:k
indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
matTS_norm1(:, iType) = mean(matAvgTS1(:, indCell_sort{iType}), 2)';
matTS_norm2(:, iType) = mean(matAvgTS2(:, indCell_sort{iType}), 2)';
end
%

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 1200 1000 340])
imagesc(zscore(matAvgTS1(:, indCelldfl))')
colormap(hot)
set(gca, 'CLim', [-2 10])
set(gca, 'YTick', find(diff(sortedIDXdfl)>0), 'XTickLabel', 20:20:120, 'TickDir', 'out')
% print(fullfile(dirFig, 'DFL1_Clustering_Tabla_5Clusters_matTS_neg1pos8'), '-r200', '-dtiff')

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 1200 1000 200])
plot(zscore(matTS_norm1), 'LineWidth', 2)
cMap_sort = hsv(k);
cMap_sort(2,:) = [206 182 49]./255;
colororder(cMap_sort)
axis tight
set(gca, 'XTickLabel', 20:20:120, 'LineWidth', 2, 'TickDir', 'out', 'box', 'off')
% print(fullfile(dirFig, 'DFL1_Clustering_Tabla_5Clusters_avgTS_zscore'), '-depsc')


% Same groups for movie 2
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 1200 1000 340])
imagesc(zscore(matAvgTS2(:, indCelldfl))')
colormap(hot)
set(gca, 'CLim', [-2 10])
set(gca, 'YTick', find(diff(sortedIDXdfl)>0), 'XTickLabel', 20:20:120, 'TickDir', 'out')
% print(fullfile(dirFig, 'DFL1_Clustering_Tabla_5Clusters_matTS_neg1pos8'), '-r200', '-dtiff')

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 1200 1000 200])
plot(zscore(matTS_norm2), 'LineWidth', 2)
cMap_sort = hsv(k);
cMap_sort(2,:) = [206 182 49]./255;
colororder(cMap_sort)
axis tight
set(gca, 'XTickLabel', 20:20:120, 'LineWidth', 2, 'TickDir', 'out', 'box', 'off')
% print(fullfile(dirFig, 'DFL1_Clustering_Tabla_5Clusters_avgTS_zscore'), '-depsc')



figure;
set(gcf, 'Color', 'w')
cMap_sort = hsv(k);
cMap_sort(2,:) = [206 182 49]./255
% cMap_sort = [215 25 28; 253 174 97; 255 255 191; 171 217 233; 44 123 182]./255;
% cMap_sort = [208 28 139; 241 182 218; 247 247 247; 184 225 134; 77 172 38]./255;
[d1 d2] = size(infoCells(1).imgFOV);
% imagesc(zeros(size(infoCells(1).imgFOV)));
% colormap(gray)
% imagesc(imgFOV);
% colormap(sp(3), gray);
% hold on;
for iType = 1:k
for iC = 1:size(indCell_sort{iType}, 1)
Coor = cellPix(indCellValid(indCell_sort{iType}(iC, 1))).contourCell{1};
plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iType, :)); hold on;
text(Coor(1,end), Coor(2,end), num2str(indCellValid(indCell_sort{iType}(iC, 1))), ...
        'color', 'k')
end
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20]) %, 'Color', 'k', 'Box', 'off')
set(gca, 'XTick', [], 'YTick', [])




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
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6','Location', 'best')
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
Coor = cellPix(indCellValid(indCell_sort{iType}(iC, 1))).contourCell{1};
plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iType, :)); hold on;
% text(Coor(1,end), Coor(2,end), num2str(indCell_sort{iType}(iC, 1)), ...
%         'color', 'k')
end
end
axis on
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20], 'Color', 'k')




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
% Coor = neuron.get_contours(thr); 
% CC = Coor;
% for i = 1:size(Aor,2)
%     %         cont = medfilt1(Coor{i}')';
%     cont = Coor{i};
%     if size(cont,2) > 1
%         plot(cont(1,1:end),cont(2,1:end),'Color',cmap(i+size(Aor,2),:), 'linewidth', ln_wd); hold on;
%     end
% end





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

