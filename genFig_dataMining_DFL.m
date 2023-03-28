% genFig_dataMining_DFL.m
%
% 2022/11/01 SHP: started new by applying longitudinally registered cells,
% instead of looking at session to session data
% 2021/12/21 SHP: working on the part that selects good cells
% 2021/09/08 SHP
% Digging DFL data to find out something
% started from "genFig_functionalMapFOV.m"

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
fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_DFL, 'cellTS', 'cellPix')
cellTS_DFL = cellTS;
clear cellTS

fname_caTSFOV_BPM = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_BPM, 'cellTS')
cellTS_BPM = cellTS;
clear cellTS

fname_caTSFOV_RS = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_RS, 'cellTS_RS', 'resultsRS')

%%
indCellValid = find(cat(1, cellTS_DFL.nTrial1_total)>8); % 
% indCellValid_trial = find(cat(1, cellTS.nTrial1_total)>8); % cells that have more than 8 trials for movie 1
% 
for ii = 1:length(cellTS_DFL)
    minsnrmovie1(ii, 1) = min(cellTS_DFL(ii).snr_movie1);
end

% indCellValid_snr = find(minsnrmovie1>0.1);   

% indCellValid = intersect(indCellValid_trial, indCellValid_snr);

clear matAvg*
for iCell = 1:length(indCellValid)   
    
    matAvgTS1(:, iCell) = mean(cellTS_DFL(indCellValid(iCell)).matTS_movie1)'; % now it's scaled dF
    matAvgTS2(:, iCell) = mean(cellTS_DFL(indCellValid(iCell)).matTS_movie2)'; %
    
    steAvgTS1(:, iCell) = std((cellTS_DFL(indCellValid(iCell)).matTS_movie1)./sqrt(size(cellTS_DFL(indCellValid(iCell)).matTS_movie1, 1)-1))'; % now it's scaled dF
    steAvgTS2(:, iCell) = std((cellTS_DFL(indCellValid(iCell)).matTS_movie2)./sqrt(size(cellTS_DFL(indCellValid(iCell)).matTS_movie2, 1)-1))'; 

end

%% DFL averaged signal corr
tsDFL_avg = struct([]);
for iCell = 1:length(cellTS_DFL)
    tsDFL_avg(iCell).avgTS1 = mean(zscore(cellTS_DFL(iCell).matTS_movie1, 0, 2))'; %
    tsDFL_avg(iCell).avgTS2 = mean(zscore(cellTS_DFL(iCell).matTS_movie2, 0, 2))'; %
    tsDFL_avg(iCell).steTS1 = std((cellTS_DFL(iCell).matTS_movie1)./sqrt(size(cellTS_DFL(iCell).matTS_movie1, 1)-1))'; 
    tsDFL_avg(iCell).steTS2 = std((cellTS_DFL(iCell).matTS_movie2)./sqrt(size(cellTS_DFL(iCell).matTS_movie2, 1)-1))'; 
end



%% Drawboard for trial-to-trial correlation computation and gather the corr values
% compute RS and DFL pairwise corr (noise+signal corr) for all the sessions
for iS = 1:length(setDateSession)
    clear matR
    [matR, matP] = corr(resultsRS(iS).C_raw', 'type', 'Spearman');
    resultsCorr(iS).matR_RS = matR;
    
    dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    load(fullfile(dirProcdata_session, 'DFL_ts_tML'));

    clear catMatR
    for iTrial = 1:size(tS_session(1).matTS_C_raw_zscore, 3)
        clear matR        
        [matR] = corr(tS_session(1).matTS_C_raw_zscore(:,:,iTrial), 'type', 'Spearman');
        
        catMatR(:,:,iTrial) = matR;
    end
    resultsCorr(iS).matR_DFL{1} = catMatR;
    resultsCorr(iS).matR_DFL_avg{1} = corr(tS_session(1).avgTS_C_raw_zscore, 'type', 'Spearman');
    
    clear catMatR
    for iTrial = 1:size(tS_session(2).matTS_C_raw_zscore, 3)
        clear matR        
        [matR] = corr(tS_session(2).matTS_C_raw_zscore(:,:,iTrial), 'type', 'Spearman');
        
        catMatR(:,:,iTrial) = matR;
    end
    resultsCorr(iS).matR_DFL{2} = catMatR;
    resultsCorr(iS).matR_DFL_avg{2} = corr(tS_session(2).avgTS_C_raw_zscore, 'type', 'Spearman');

end

% For all the longitudinally tracted cells
indCellValid = find(cat(1, cellTS_DFL.nTrial1_total)>8); % 
cellIDAcrossDay_validCell = cat(2, cellIDAcrossDay(indCellValid, :), indCellValid);

setPair = nchoosek(1:length(indCellValid), 2);
% later I can reorganize this into matrix using the below line
% ind = sub2ind([length(indCellValid) length(indCellValid)], setPair(:,1), setPair(:, 2));
registeredCellPairCorr = struct([]);

for iPair = 1:length(setPair)
    
    idPair = setPair(iPair,:);
    
    setSession = find(sum(~isnan(cellIDAcrossDay_validCell(setPair(iPair,:), 1:end-1)))==2);
    
       
    corrRS = []; corrDFL_1 = []; corrDFL_2 = []; corrDFL_1_avgSes = []; corrDFL_2_avgSes = [];
    for iSes = 1:length(setSession) %size(cellIDAcrossDay_validCell, 2)-1
%         if sum(isnan(cellIDAcrossDay_validCell(setPair(iPair,:), iSes))) > 0
%             continue;
%         end
        idSes = setSession(iSes);

        indCell_curPair = cellIDAcrossDay_validCell(setPair(iPair,:), idSes);
        
        corrRS(iSes,:) = [idSes resultsCorr(idSes).matR_RS(indCell_curPair(1), indCell_curPair(2))];
        corrDFL_1 = cat(1, corrDFL_1, cat(2, repmat(idSes, size(resultsCorr(idSes).matR_DFL{1}, 3), 1), ...
            squeeze(resultsCorr(idSes).matR_DFL{1}(indCell_curPair(1), indCell_curPair(2), :))));
        corrDFL_2 = cat(1, corrDFL_2, cat(2, repmat(idSes, size(resultsCorr(idSes).matR_DFL{2}, 3), 1), ...
            squeeze(resultsCorr(idSes).matR_DFL{2}(indCell_curPair(1), indCell_curPair(2), :))));
        corrDFL_1_avgSes(iSes,:) = [idSes resultsCorr(idSes).matR_DFL_avg{1}(indCell_curPair(1), indCell_curPair(2))];
        corrDFL_2_avgSes(iSes,:) = [idSes resultsCorr(idSes).matR_DFL_avg{2}(indCell_curPair(1), indCell_curPair(2))];
    end
    
    registeredCellPairCorr(iPair).idPair = idPair;
    registeredCellPairCorr(iPair).setSession = setSession;
    registeredCellPairCorr(iPair).corrRS = corrRS;
%     registeredCellPairCorr(iPair).corrRS_mean = mean(corrRS(:,2));
    registeredCellPairCorr(iPair).corrDFL_1 = corrDFL_1;
    registeredCellPairCorr(iPair).corrDFL_2 = corrDFL_2;
    registeredCellPairCorr(iPair).corrDFL_1_avgSes = corrDFL_1_avgSes;
    registeredCellPairCorr(iPair).corrDFL_2_avgSes = corrDFL_2_avgSes;
    
    %% movie signal correlation using averaged time series
    R1 = corr(tsDFL_avg(indCellValid(idPair(1))).avgTS1, tsDFL_avg(indCellValid(idPair(2))).avgTS1, 'type', 'Spearman');
    registeredCellPairCorr(iPair).corrDFL_1_signal = R1;
    R2 = corr(tsDFL_avg(indCellValid(idPair(1))).avgTS2, tsDFL_avg(indCellValid(idPair(2))).avgTS2, 'type', 'Spearman');
    registeredCellPairCorr(iPair).corrDFL_2_signal = R2;
    
end
        
for iPair = 1:length(registeredCellPairCorr)
    
    if isempty(registeredCellPairCorr(iPair).setSession)
        setCorrMean(iPair, 1:5) = 0; %NaN;
        continue;
    end
    
    setCorrMean(iPair, 1) = mean(registeredCellPairCorr(iPair).corrRS(:,2));
    setCorrMean(iPair, 2) = mean(registeredCellPairCorr(iPair).corrDFL_1(:,2));
    setCorrMean(iPair, 3) = mean(registeredCellPairCorr(iPair).corrDFL_2(:,2));
    setCorrMean(iPair, 4) = mean(registeredCellPairCorr(iPair).corrDFL_1_avgSes(:,2));
    setCorrMean(iPair, 5) = mean(registeredCellPairCorr(iPair).corrDFL_2_avgSes(:,2));
    setCorrMean(iPair, 6) = mean(registeredCellPairCorr(iPair).corrDFL_1_signal);
    setCorrMean(iPair, 7) = mean(registeredCellPairCorr(iPair).corrDFL_2_signal);
    
end
        

ind = sub2ind([length(indCellValid) length(indCellValid)], setPair(:,1), setPair(:, 2));
matCorrRS = NaN(length(indCellValid));
matCorrRS(ind) = setCorrMean(:,1);
matCorrMov1 = NaN(length(indCellValid));
matCorrMov1(ind) = setCorrMean(:,2);
matCorrMov2 = NaN(length(indCellValid));
matCorrMov2(ind) = setCorrMean(:,3);
matCorrMov1_signal = NaN(length(indCellValid));
matCorrMov1_signal(ind) = setCorrMean(:,6);
matCorrMov2_signal = NaN(length(indCellValid));
matCorrMov2_signal(ind) = setCorrMean(:,7);

tempSubset = 1:109;
figure;
set(gcf, 'Color', 'w')
sp(1) = subplot(1,3,1);
imagesc(matCorrMov1(tempSubset, tempSubset))
title('mov 1')
sp(2) = subplot(1,3,2);
imagesc(matCorrMov1_signal(tempSubset, tempSubset))
title('mov 1:signal')
sp(3) = subplot(1,3,3);
imagesc(matCorrRS(tempSubset, tempSubset))
title('RS')
set(sp, 'CLim', [-1 1].*0.6)
colormap('jet')


%%
[matR_sort, sortedRow] = sort(setCorrMean, 'descend'); % 1st column: RS, 2nd column: DFL1
% setBothHigh = intersect(sortedRow(1:100,1), sortedRow(1:100,2));
% setPair(setBothHigh(1:20),:)

% top 20 pairs of highest postivie movie 1 correlation
setMV = indCellValid(setPair(sortedRow(1:20,6),:));
fig_MVcorrpairs = figure;
set(fig_MVcorrpairs, 'Color', 'w', 'Position', [150 150 570 413])
dim_fov = size(infoCells(1).imgFOV);
for iP = 1:length(setMV)
    
    curP = setMV(iP,:);
    curP_centercoords = cat(1, cellPix(curP(1)).centerCell(1,:), cellPix(curP(2)).centerCell(1,:)); % in image coords (row, column)
        
    figure(fig_MVcorrpairs);
    hold on;
    plot(curP_centercoords(1,2), curP_centercoords(1,1), 'b.', 'MarkerSize', 10); hold on;  
    text(curP_centercoords(1,2)+1, curP_centercoords(1,1), num2str(curP(1)));  
    plot(curP_centercoords(2,2), curP_centercoords(2,1), 'b.', 'MarkerSize', 10); hold on;
    text(curP_centercoords(2,2)+1, curP_centercoords(2,1), num2str(curP(2)));
    line(curP_centercoords(:,2), curP_centercoords(:,1), 'Color', 'b');
    set(gca, 'YDir', 'reverse', 'XLim', [0-10 dim_fov(2)+10], 'YLim', [0-10 dim_fov(1)+10], 'TickDir', 'out');
end
title(sprintf('%s: top 20 pairs with highest positive correlation during movie 1', nameSubj))


figure
subplot(3,1,1)
imagesc(zscore(cellTS_DFL(193).matTS_movie1')')
subplot(3,1,2)
imagesc(zscore(cellTS_DFL(124).matTS_movie1')')
subplot(3,1,3)
imagesc(zscore(cellTS_DFL(63).matTS_movie1')')


%% Residuals
% for each time bin (e.g. 1s), gather the avg activity of cell 1 and
% correlate them with cell 2, so pair-wise correlation for each time bin
% then see whether there are some correlation change across time

% For all the longitudinally tracted cells
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 
cellIDAcrossDay_validCell = cat(2, cellIDAcrossDay(indCellValid, :), indCellValid);

setPair = nchoosek(1:length(indCellValid), 2);

sizeBin_msec = 1000; %10;
sizeBin = sizeBin_msec./100; %100ms sampling (10hz imaging)

for iPair = 1:length(setPair)
    
    idPair = setPair(iPair,:); % this is an index of "cells more than 8 trials" (i.e. 1 - 157 for tabla, 1-7 for max)
    idPair_cellTS = cellIDAcrossDay_validCell(idPair, end); % this is the index of "cellTS" (i.e. all the cells longitudinally registered)
    
    setSession = find(sum(~isnan(cellIDAcrossDay_validCell(idPair, 1:end-1)))==2);
    
    % cell 1 of this pair
    cellTS(idPair_cellTS(1)).idAcrossSession
    cellTS(idPair_cellTS(1)).nTrial1_set
    
    temp = cat(1, 1, cellTS(idPair_cellTS(1)).nTrial1_set(1:end-1)+1);
    setTrial1 = cat(2, temp, cellTS(idPair_cellTS(1)).nTrial1_set);
    
    locSes = ismember(cellTS(idPair_cellTS(1)).idAcrossSession(:,1), setSession);
    
    % cell 2 of this pair
    
    
%     for iSes = 1:length(setSession) %
%         idSes = setSession(iSes);
%         
%         
        

    
    iCell = 1;
tMat = cellTS(indCellValid(iCell)).matTS_movie1(:, 1:1200)';

iBinT = 1;
(iBinT-1)*sizeBin+1:iBinT*sizeBin
    
    
       
    corrRS = []; corrDFL_1 = []; corrDFL_2 = []; corrDFL_1_avgSes = []; corrDFL_2_avgSes = [];
    for iSes = 1:length(setSession) %size(cellIDAcrossDay_validCell, 2)-1
%         if sum(isnan(cellIDAcrossDay_validCell(setPair(iPair,:), iSes))) > 0
%             continue;
%         end
        idSes = setSession(iSes);

        indCell_curPair = cellIDAcrossDay_validCell(setPair(iPair,:), idSes);
        
        corrRS(iSes,:) = [idSes resultsCorr(idSes).matR_RS(indCell_curPair(1), indCell_curPair(2))];
        corrDFL_1 = cat(1, corrDFL_1, cat(2, repmat(idSes, size(resultsCorr(idSes).matR_DFL{1}, 3), 1), ...
            squeeze(resultsCorr(idSes).matR_DFL{1}(indCell_curPair(1), indCell_curPair(2), :))));
        corrDFL_2 = cat(1, corrDFL_2, cat(2, repmat(idSes, size(resultsCorr(idSes).matR_DFL{2}, 3), 1), ...
            squeeze(resultsCorr(idSes).matR_DFL{2}(indCell_curPair(1), indCell_curPair(2), :))));
        corrDFL_1_avgSes(iSes,:) = [idSes resultsCorr(idSes).matR_DFL_avg{1}(indCell_curPair(1), indCell_curPair(2))];
        corrDFL_2_avgSes(iSes,:) = [idSes resultsCorr(idSes).matR_DFL_avg{2}(indCell_curPair(1), indCell_curPair(2))];
    end
    
    registeredCellPairCorr(iPair).idPair = idPair;
    registeredCellPairCorr(iPair).setSession = setSession;
    registeredCellPairCorr(iPair).corrRS = corrRS;
%     registeredCellPairCorr(iPair).corrRS_mean = mean(corrRS(:,2));
    registeredCellPairCorr(iPair).corrDFL_1 = corrDFL_1;
    registeredCellPairCorr(iPair).corrDFL_2 = corrDFL_2;
    registeredCellPairCorr(iPair).corrDFL_1_avgSes = corrDFL_1_avgSes;
    registeredCellPairCorr(iPair).corrDFL_2_avgSes = corrDFL_2_avgSes;
end


% compute residual & mean
for iCell = 1:length(indCellValid)

tMat = cellTS(indCellValid(iCell)).matTS_movie1;
tMatNorm = zscore(cellTS(indCellValid(iCell)).matTS_movie1')';
tMatNorm_mean = mean(tMatNorm);
tMatNorm_res = tMatNorm - tMatNorm_mean;

figure(100);clf;
subplot(2,1,1)
imagesc(tMatNorm)
title(sprintf('Cell #%d, ID%d: before & after mean subtraction', iCell, indCellValid(iCell)))
set(gca, 'XTickLabel', 20:20:120)
subplot(2,1,2)
imagesc(tMatNorm_res)
set(gca, 'CLim', [0 10])
set(gca, 'XTickLabel', 20:20:120)

figure(200);clf;
subplot(2,1,1)
plot(tMatNorm(1:5, :)')
hold on;
plot(tMatNorm_mean, 'k-', 'LineWidth', 2)
axis tight
title(sprintf('Cell #%d, ID%d: before & after mean subtraction (first 5 trials)', iCell, indCellValid(iCell)))
set(gca, 'XTickLabel', 20:20:120)
subplot(2,1,2)
plot(tMatNorm_res(1:5, :)')
axis tight
set(gca, 'XTickLabel', 20:20:120)

input('')
end


    
%     iC = setPair(iPair, 1);
%     jC = setPair(iPair, 2);
    
  

iCell = 1; % size(cellIDAcrossDay_validCell, 1)
setSession = find(~isnan(cellIDAcrossDay_validCell(iCell, 1:size(cellIDAcrossDay_validCell, 2)-1)));
iSession = 1;
idSession = setSession(iSession);

% indices of valid members of cells from this session
indCellSession = cellIDAcrossDay_validCell(~isnan(cellIDAcrossDay_validCell(:, idSession)), idSession);


% pairwise correlation during movie

%     resultsCorr_DFL(iS).matR = matR;




% end


%% PCA part - now separated and saved into "analCa_DFL_PCA.m"
% [coeff, score, latent, tsquared, explained] = pca(zscore(matAvgTS1)');
% [sortedScore, indCell] = sort(score(:,1:10), 'descend');
% 
% 
% % [sortedScore, indCell] = sort(score(:,1), 'descend');
% % [sortedScore2, indCell2] = sort(score(:,2), 'descend');
% 
% explained(1:10)
% 
% figure;
% subplot(2,1,1)
% imagesc(zscore(matAvgTS1)');
% set(gca, 'CLim', [0 10])
% colormap(hot)
% title('movie 1 responses (zscore)')
% subplot(2,1,2)
% imagesc(zscore(matAvgTS1(:,indCell(:,1)))');
% set(gca, 'CLim', [0 10])
% colormap(hot)
% title('movie 1 responses (zscore): sorted along PC1')
% 
% figure;
% plot(coeff(:,1:3))
% legend('PC1', 'PC2', 'PC3')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% in the middle of working
% 
% [sortedScore, indCell] = sort(score(:,1:10), 'descend');
% 
% resultsPCA_block(iB, iMovie).explained = explained;
% resultsPCA_block(iB, iMovie).coeff = coeff(:, 1:10);
% resultsPCA_block(iB, iMovie).score = score(:, 1:10);
% resultsPCA_block(iB, iMovie).indCellSorted = indCell;
% 
% resultsPCA(iSubj).resultsPCA_block = resultsPCA_block;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
% 
% 
% 
% [sortedScore_abs, indCell_abs] = sort(abs(score(:,1)), 'descend');
% 
% figure;
% % cMap_sort = jet(length(indCell)); %jet(length(indCell)); %hsv(k);
% cMap_sort = autumn(length(indCell));
% 
% [d1 d2] = size(infoCells(1).imgFOV);
% % imagesc(imgFOV); 
% % colormap(sp(3), gray);
% % hold on;
% for iCell = 1:length(indCell)
%     iC = indCell(iCell);
%         Coor = cellPix(indCellValid(iC)).contourCell{1};
%         plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iCell, :)); hold on;
% end
% set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
% 
% 
% k = 4;
% [IDXdfl, C, SUMD] = kmeans(zscore(matAvgTS1)', k, 'Distance', 'correlation');
% [sortedIDXdfl, indCelldfl] = sort(IDXdfl);
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
% end
% 
% % 
% % Plotting
% fig_summary_DFL = figure;
% set(gcf,  'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1085 750])
% clear sp
% % 
% % % 1. averaged amplitude for each condition
% % figure(fig_summary_DFL);
% % sp(1) = subplot('Position', [0.2 0.65 0.75 0.3]);
% % imagesc(catAvgMatTS(:, indCelldfl)');
% % colormap(hot);
% % set(sp(1), 'CLim', [-1 1].*4)
% % set(sp(1), 'YTick', find(diff(sortedIDXdfl)>0), 'YTickLabel', [])
% % set(sp(1), 'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
% % set(sp(1), 'TickDir', 'out')
% % box off
% % colorbar;
% % title(sprintf('%s %s: averaged response to movie %s', nameSubj, dateSession, tS_session(iMovie).idStim))
% % ylabel('Cells (sorted)')
% % xlabel('Time (s)')
% % 
% % 2. Clustering results on 2-d PC space
% figure(fig_summary_DFL);
% sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% for iType = 1:k
%         plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
%         hold on;
% end
% legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Location', 'best')
% xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
% ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
% set(gca, 'TickDir', 'out')
% box off
% axis square
% title('Clustering based on movie response on PC space')
% 
% % 3. Clustering results on imaging field of view
% figure(fig_summary_DFL);
% sp(3) = subplot('Position', [0.55 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% [d1 d2] = size(infoCells(1).imgFOV);
% % imagesc(imgFOV); 
% % colormap(sp(3), gray);
% % hold on;
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         Coor = cellPix(indCellValid(indCell_sort{iType}(iC, 1))).contourCell{1};
%         plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iType, :)); hold on;
%     end
% end
% axis on
% set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
% 
% 
% 
% %% PCAs on 20-s blocks
% % block information
% paramCond.condName = {'Scene Motion', 'Human Faces', 'Monkey Bodies', 'Monkey Faces', 'Object Motion', 'Optic Flow'}';
% paramCond.condOrder = [1 2 6 3 5 4; 5 6 4 1 3 2]';
% [b, indCondReorder] = sort(paramCond.condOrder, 1);
% 
% indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 
% % indCellValid_trial = find(cat(1, cellTS.nTrial1_total)>8); % cells that have more than 8 trials for movie 1
% % 
% for ii = 1:length(cellTS)
%     minsnrmovie1(ii, 1) = min(cellTS(ii).snr_movie1);
% end
% 
% % indCellValid_snr = find(minsnrmovie1>0.1);   
% 
% % indCellValid = intersect(indCellValid_trial, indCellValid_snr);
% 
% clear matAvg*
% for iCell = 1:length(indCellValid)   
%     
%     matAvgTS1(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie1)'; % now it's scaled dF
%     matAvgTS2(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie2)'; %
%     
%     steAvgTS1(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie1)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie1, 1)-1))'; % now it's scaled dF
%     steAvgTS2(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie2)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie2, 1)-1))'; 
% 
% end
% 
% matAvgTS1_block = reshape(matAvgTS1(1:1200, :)', size(matAvgTS1, 2), 200, 6);
% matAvgTS2_block = reshape(matAvgTS2(1:1200, :)', size(matAvgTS2, 2), 200, 6);
% 
% matAvgTS_block_reorder{1} = matAvgTS1_block(:,:,indCondReorder(:,1));
% matAvgTS_block_reorder{2} = matAvgTS2_block(:,:,indCondReorder(:,2));
% % matAvgTS2_block_reorder = matAvgTS2_block(:,:,indCondReorder(:,2));
% 
% figure
% for iR = 1:6
% subplot(2, 6, iR);
% imagesc(zscore(matAvgTS2_block(:,:,iR)')')
% end
% for iR = 1:6
% subplot(2,6,iR+6);
% imagesc(zscore(matAvgTS_block_reorder{2}(:,:,iR)')')
% end
% 
% for iMovie = 1:2
%     for iB = 1:6
%         
%         [coeff, score, latent, tsquared, explained] = pca(zscore(matAvgTS_block_reorder{iMovie}(:,:,iB)')'); %pca(zscore(matAvgTS1)');
%         [sortedScore, indCell] = sort(score(:,1:10), 'descend');
%         
%         %     explained(1:10)
%         
% %             figure;
% %             set(gcf, 'Color', 'w', 'Position', [700 700 885 415])
% %             subplot(1,3,1);
% %             plot(coeff(:,1:3))
% %             legend('PC1', 'PC2', 'PC3')
% %             set(gca, 'XTick', 0:100:200, 'XTickLabel', 0:10:20)
% %             xlabel('Time (s)')
% %             title(sprintf('%s', paramCond.condName{iB}))
% %             subplot(1,3,2)
% %             imagesc(zscore(matAvgTS_block_reorder{iMovie}(indCell(:,1),:,iB)')')
% %             colormap(hot)
% %             set(gca, 'XTick', 0:50:200, 'XTicklabel', 0:5:20)
% %             title('PC1 sorted')
% %             xlabel('Time (s)')
% %             ylabel('Cells (sorted)')
% %             subplot(1,3,3)
% %             imagesc(zscore(matAvgTS_block_reorder{iMovie}(indCell(:,2),:,iB)')')
% %             colormap(hot)
% %             set(gca, 'XTick', 0:50:200, 'XTicklabel', 0:5:20)
% %             title('PC2 sorted')
% %             xlabel('Time (s)')
% %             ylabel('Cells (sorted)')
% % 
% %             print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_DFLmovie%d_PCA_BlockID%d', nameSubj, FOV_ID, iMovie, iB)), '-depsc')
%         
%         resultsPCA_block(iB, iMovie).explained = explained;
%         resultsPCA_block(iB, iMovie).coeff = coeff(:, 1:10);
%         resultsPCA_block(iB, iMovie).score = score(:, 1:10);
%         resultsPCA_block(iB, iMovie).indCellSorted = indCell;
%         
%     end
% end
% 
% resultsPCA(iSubj).resultsPCA_block = resultsPCA_block;
% 
% %%%% next step is to draw FOV pc score
% for iB = 1:size(resultsPCA_block, 1)
%     fig_fov = figure;
%     set(fig_fov, 'Color', 'w', 'Position', [1200 560 1060 390])
%     for iMovie = 1:2
%         iPC = 1;
%         curIndCell = resultsPCA_block(iB, iMovie).indCellSorted(:,iPC);
%         cMap_sort = jet(length(curIndCell)); %cool(length(curIndCell)); %jet(length(indCell)); %hsv(k);
%         setSortedIndCellValid = indCellValid(curIndCell);
%         
%         figure(fig_fov);
%         subplot(1,2,iMovie);
%         for iC = 1:length(setSortedIndCellValid)
%             
%             idCurCell = setSortedIndCellValid(iC);
%             plot(cellPix(idCurCell).centerCell(1,2)+shifts(cellPix(idCurCell).idAcrossSession(1,1),1), ...
%                 cellPix(idCurCell).centerCell(1,1)+shifts(cellPix(idCurCell).idAcrossSession(1,1),2), 'o',...
%                 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap_sort(iC, :)); hold on;
%             
%             %         Coor = cellPix(indCellValid(iC)).contourCell{1};
%             %         plot(Coor(1,:), Coor(2,:), '.', 'Color', cMap_sort(iCell, :)); hold on;
%         end
%         [d1, d2] = size(infoCells(1).imgFOV);
%         set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20])
%         title(sprintf('%s: movie%d', paramCond.condName{iB}, iMovie))
%         
%     end
% end
% 
% matColorIndex = max(min(fix(((curMatR - cmin)./(cmax-cmin)*256)+1), 256), 1); % saturating it to cmax and cmin
% 
% 
%     figMap_RS = figure;
%     set(figMap_RS, 'Color', 'w')
% for iCell = 1:length(validCell)
%     
%     idCurCell = validCell(iCell);
%     lw = 2;
%     if idCurCell == curCells_id(iSetCell)
%         lw = 4;
%     end
%     
%     figure(figMap_RS);
%     subplot(nRow, nCol, iSetCell);
%     plot(cellPix(idCurCell).centerCell(1,2)+shifts(curCells_session(iSetCell),1), ...
%                 infoCells(curCells_session(iSetCell)).cellCenter(idCurCell,1)++shifts(curCells_session(iSetCell),2), 'o',...
%                 'MarkerSize', 8, 'LineWidth', lw, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap_bwr(matColorIndex(idCurCell), :));
%     hold on;
% end
% [d1 d2] = size(infoCells(1).imgFOV);
% axis equal
% set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20]);
% title(sprintf('Resting State corr: %s session%d Cell%d', nameSubj, curCells_session(iSetCell), curCells_id(iSetCell)))
% 
% 
% 
% % for iS = 1:length(setDateSession)
% %     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
% %     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% %     load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
% %     
% % %     resultsDFL(iS).tS_session = tS_session;
% %     
% %     [matR, matP] = corr(tS_session(1).avgTS_C_raw_zscore, 'type', 'Spearman');
% %     resultsCorr_DFL(iS).matR = matR;
% % end


%% spatial map
tempA = cat(2, cellPix(indCellValid).repPix);
tempA(~isnan(tempA)) = 1;
% tempA(:, indCellValid) = tempA(:, indCellValid).*10;

imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));

figure;
set(gcf, 'Color', 'w')
imagesc(imgCells_2d)
colormap(turbo)


% pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);


%% more serious clustering attempt
[a, c, totalSS] = kmeans(zscore(matAvgTS1)', 1, 'Distance', 'correlation'); % cell
Clustering.totalSS_cell = totalSS;
[a, c, totalSS] = kmeans(zscore(matAvgTS1), 1, 'Distance', 'correlation'); % time
Clustering.totalSS_time = totalSS;

setK = 2:10;
numReplicates = 5;
opts = statset('Display','final');
numRepeat = 100;
for iK = 1:length(setK)
    
    K = setK(iK);                
    
    % max of abs
    cell_indCluster = NaN(size(matAvgTS1, 2), numRepeat);
    cell_sumD = NaN(K, numRepeat);
    
    time_indCluster = NaN(size(matAvgTS1, 1), numRepeat);
    time_sumD = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; maxabs ::\n', K);
        
        % Cluster single units based on whole brain correlation
        [IDX_SU, C_SU, SUMD_SU] = kmeans(zscore(matAvgTS1)', K,...
            'Replicates', numReplicates, 'Options', opts, 'Distance', 'correlation');
        
        cell_indCluster(:, iRep) = IDX_SU;
        cell_sumD(:, iRep) = SUMD_SU;
        
        % Cluster ROIs based on single unit correlation
        [IDX_time, C_time, SUMD_time] = kmeans(zscore(matAvgTS1), K,...
            'Replicates', numReplicates, 'Options', opts, 'Distance', 'correlation');
        
        time_indCluster(:, iRep) = IDX_time;
        time_sumD(:, iRep) = SUMD_time;
        
    end

    Clustering.resultKMeans(iK).cell_indCluster = cell_indCluster;
    Clustering.resultKMeans(iK).cell_sumD = cell_sumD;
    Clustering.resultKMeans(iK).time_indCluster = time_indCluster;
    Clustering.resultKMeans(iK).time_sumD = time_sumD;

end

paramClustering.methods = 'KMeans using Correlation distance';
paramClustering.setK = setK;
paramClustering.numReplicates = numReplicates;
paramClustering.numRepeat = numRepeat;


%% Explained variance elbow plot
setK = paramClustering.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
%     indClust = Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
%     [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(:,iK) = sum(Clustering.resultKMeans(iK).cell_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering.totalSS_cell;
% [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
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
    plot(iK+1, mean(Clustering.resultKMeans(iK).cell_sumD, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')

%%
% clear all;
% 
% ss = pwd;
% if ~isempty(strfind(ss, 'Volume')) % if it's local
%     dirProjects = '/Volumes/NIFVAULT/projects/parksh';
%     dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
%     dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
% else % on virtual machine
%     dirProjects = '/nifvault/projects/parksh';
%     dirProcdata = '/nifvault/procdata/parksh';
%     dirRawdata = '/nifvault/rawdata/parksh';
% end
% 
% % setNameSubj = {'Tabla', 'Max'};
% flagSavePPTX = 0; %1;
% 
% % get session info
% nameSubj = 'Tabla';
% FOV_ID = 1;
% [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);
% 
% for iSession = 1:nSession
% % iSession = 1; 
% dateSession = setDateSession{iSession}; %'20191113'; %setDateSession{iSession};
% 
% dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% dirPreproc = fullfile(dirProcdata_session, '_preproc');
% 
% dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
% 
% 
% %% Read source data
% addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
% cnmfe_setup;
% d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
% 
% load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
% 
% validIndCell = [];
% validIndCell(:,1) = 1:length(neuron.ids);
% if strcmpi(nameSubj, 'max')
%     load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
%     validIndCell = indCell.validCell;
% end
% 
% 
% % get the contours and image field of view
% % neuron_b = neuron.batches{1}.neuron;
% 
% % Generate cell location map within FOV
% thr = 0.5; % the lower the smaller (more centralized) the contour
% cellColor = [1 1 1];
% widthContour = 1;
% [d1,d2] = size(neuron.Cn);
% 
% figure;
% imagesc(zeros(d1, d2)); % background
% colormap(gray);
% caxis([0 0.1]);
% hold on;
% 
% CC = cell(size(neuron.A, 2),1);
% CR = cell(size(neuron.A, 2),2);
% for i = 1:size(neuron.A ,2)
%     A_temp = full(reshape(neuron.A(:,i),d1,d2));
%     A_temp = medfilt2(A_temp,[3,3]);
%     A_temp = A_temp(:);
%     [temp,ind] = sort(A_temp(:).^2,'ascend');
%     temp =  cumsum(temp);
%     ff = find(temp > (1-thr)*temp(end),1,'first');
%     if ~isempty(ff)
%         CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
%         fp = find(A_temp >= A_temp(ind(ff)));
%         [ii,jj] = ind2sub([d1,d2],fp);
%         CR{i,1} = [ii,jj]';
%         CR{i,2} = A_temp(fp)';
%     end
%     hold on;
% end
% axis off
% title(sprintf('%s: %s', nameSubj, dateSession))
% % %save
% % print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');
% 
% % end
% %
% Coor = neuron.get_contours(thr); % Coor = get_contours(obj, thr, ind_show); % ind_show: indices of cells you want to get contours
% imgFOV = neuron.Cn.*neuron.PNR;
% 
% % % draw all the contours
% % neuron_b.show_contours([], [], imgFOV, 'true');
% 
% figure
% [center] = neuron.estCenter();
% center = center(validIndCell, :);
% imagesc(imgFOV)
% hold on
% plot(center(:,2), center(:, 1), 'r.')
% text(center(:,2)+1, center(:,1), num2str([1:size(center,1)]'), 'Color', 'w');
% axis off
% 
% 
% %% Cell quality check using SNR and PNR
% snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);
% pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2);
% 
% snrs = snrs(validIndCell);
% pnrs = pnrs(validIndCell);
% 
% [a, sortedID_snr] = sort(snrs, 'descend');
% [aa, sortedID_pnr] = sort(pnrs, 'descend');
% 
% % SNR distribution 
% figure
% histogram(snrs, 30)
% median(snrs)
% set(0, 'defaultfigurecolor', [1 1 1])
% xlabel('SNR')
% title(sprintf('%s: %s', nameSubj, dateSession))
% line([median(snrs) median(snrs)], get(gca, 'ylim'), 'Color', 'r')
% 
% % Location of a particular set cells
% setCells_rank = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr); %1:10
% figure;
% imagesc(imgFOV);
% colormap(gray);
% hold on
% plot(center(sortedID_snr(setCells_rank), 2), center(sortedID_snr(setCells_rank), 1), 'r.')
% hold on
% text(center(sortedID_snr(setCells_rank), 2)+1, center(sortedID_snr(setCells_rank), 1)+1, num2str([setCells_rank]'), 'Color', 'r');
% text(center(sortedID_pnr(setCells_rank), 2)+3, center(sortedID_pnr(setCells_rank), 1)+3, num2str([setCells_rank]'), 'Color', 'g');
% plot(center(sortedID_pnr(setCells_rank), 2), center(sortedID_pnr(setCells_rank), 1), 'g.')
% axis off
% 
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(100);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_snr(iCell)), 1:nTime)+10*(iCell-1), '-')
% tY(iCell) = mean(neuron.C_raw(validIndCell(sortedID_snr(iCell)), 1:nTime)+10*(iCell-1));
% end
% title('SNR-based sorting')
% set(gca, 'YTick', tY, 'YTickLabel', 1:nCell, 'TickDir', 'out')
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% xlabel('Time (s)')
% ylabel('Cell order')
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(200);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+10*(iCell-1), '-')
% tY(iCell) = mean(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+10*(iCell-1));
% end
% title('PNR-based sorting')
% set(gca, 'YTick', tY, 'YTickLabel', 1:nCell, 'TickDir', 'out')
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% xlabel('Time (s)')
% ylabel('Cell order')
% 
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(200);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+3*(iCell-1), '-', 'LineWidth', 2)
% end
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% 
% 
% 
% %% Movie-driven signal
% load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
% 
% % snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);
% % pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2);
% % 
% % [a, sortedID_snr] = sort(snrs, 'descend');
% % [aa, sortedID_pnr] = sort(pnrs, 'descend');
% 
% setCell = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr); %1:10; %11:20; %10; 
% tY = []; ttY = [];
% figMovie_snr = figure;
% set(figMovie_snr, 'Position', [675    31   660   930], 'name', sprintf('%s %s: PNR based sorting', nameSubj, dateSession))
% SP(1) = subplot(1,2,1);
% SP(2) = subplot(1,2,2);
% hold(SP(:), 'on');
% 
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell);
%     figure(figMovie_snr);
%     plot(SP(1), squeeze(tS_session(1).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1), '-'); 
%     tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1)));
%     
%     plot(SP(2), squeeze(tS_session(2).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1), '-'); 
%     ttY(iCell) = mean(mean(squeeze(tS_session(2).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1)));
% end
% set(SP, 'YTick', tY, 'YTickLabel', setCell, 'TickDir', 'out')
% set(SP, 'XTick', 200:200:1200,  'XTickLabel', [200:200:1200]./10);
% axis(SP, 'tight')
% xlabel(SP, 'Time (s)')
% ylabel(SP(1), 'Cell')
% title(SP(1), 'Movie 1')
% title(SP(2), 'Movie 2')
% 
% %pnr
% setCell = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr);  %1:10; %11:20; %10; 
% tY = []; ttY = [];
% figMovie_pnr = figure;
% set(figMovie_pnr, 'Position', [675    31   660   930], 'name', sprintf('%s %s: PNR based sorting', nameSubj, dateSession))
% SP(1) = subplot(1,2,1);
% SP(2) = subplot(1,2,2);
% hold(SP(:), 'on');
% 
% scalefac = 10;
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell);
%     figure(figMovie_pnr);
%     plot(SP(1), squeeze(tS_session(1).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1), '-'); 
%     tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1)));
%     
%     plot(SP(2), squeeze(tS_session(2).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1), '-'); 
%     ttY(iCell) = mean(mean(squeeze(tS_session(2).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1)));
% end
% set(SP, 'YTick', tY, 'YTickLabel', setCell, 'TickDir', 'out')
% set(SP, 'XTick', 200:200:1200,  'XTickLabel', [200:200:1200]./10);
% axis(SP, 'tight')
% xlabel(SP, 'Time (s)')
% ylabel(SP(1), 'Cell')
% title(SP(1), 'Movie 1')
% title(SP(2), 'Movie 2')
% 
% 
% %% pupil size change
% % part of eye data is not great. decided to focus on the later half of the
% % first movie, which contains body motion, object motion, face
% 
%     % get session info
%     [infoSession, opts] = readInfoSession(nameSubj);
%     S = table2struct(infoSession);
%     
%     % setExpName = {S.ExpName}';
%     setMLFilename = {S.MLFilename}';
%     
%     indDFLRuns = contains(setMLFilename, 'DFL') & cat(1, S.flagPreproc) > 0 & contains({S.stimulus}', 'set1_1'); %% containing "DFL" in filename AND flagPreproc value of 1
%     setFilename = setMLFilename(indDFLRuns);
%     
%     % Tabla
%     % 191113: eye signal not good (another range of x/y/pupil present during the first half)
%     % 191114 & 191118: eye signal good. Length of data points are 6-7ms shorter than 2min
%     % except the last run of 191118
%     % 191119: starting to lose x signal + y & pupil is bad in general
%     % 191120: x is lost, y& pupil was good in the first run, but not in the
%     % others
%     % 191121: x is lost, length of data points are 6-7ms shorter than 2min, y
%     % & pupil bad in some runs
%    % 191125: x is lost, for the first file signals look okay & length is fine 
%    
%     
%     for iFile = 1:length(setFilename)
%         filename = strcat(setFilename{iFile}, '.bhv2');
%         
%         dateSession = filename(1:6);
%         
%         if str2num(dateSession) < 191121
%             dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
%         else
%             dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
%         end
%         
%         % filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have
%         
%         %% Read the file
%         data = mlread(fullfile(dirBHV, filename)); % mlread(filename);
%         
%         %% eye data during stimulus on 
%         % Event Code Numbers & Names : TASK_START = 10; FP_ON = 20;
%         % WAIT_FOR_TR = 30; MOVIE_ON = 40; REWARD = 90; 
%         % TRIG onset = 900; TRIG offset = 990; (TTL from ML to Inscopix DAQ On & Off)
%         locStimOn = find(data.BehavioralCodes.CodeNumbers == 40);
%         time_stimOn = floor(data.BehavioralCodes.CodeTimes(locStimOn))
% %         data.AnalogData
%         
%         figure;
%         plot(data.AnalogData.Eye, '.')
%         line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
%         title(sprintf('%s: x & y gaze', filename))
%         axis tight
%         xlabel('Time')
%         ylabel('Voltage')
%         legend('X', 'Y')
%         
%         figure;
%         plot(data.AnalogData.General.Gen1, 'g.')
%         line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
%         title(sprintf('%s: pupil size change', filename))
%         axis tight
%         xlabel('Time')
%         ylabel('Voltage')
%         
%         input('')
%     end
% tempP = data.AnalogData.General.Gen1((time_stimOn:end));
% figure
% plot(tempP)
% plot(tempP(60001:120000))
% for iCell = 1:length(setCell)
% figure(200); hold on;
% idCell = setCell(iCell);
% plot(squeeze(tS_session(1).matTS(:, sortedID_pnr(idCell), 1))+scalefac*(iCell-1), '-');
% tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, sortedID_pnr(idCell), 1))+scalefac*(iCell-1)));
% end
% axis tight
% xlim([600 1200])


%% OLD CODE: 1) consistency across trials, 2) clustering of averaged movie-driven responses
% ts = struct([]);
% iMovie = 1;
% for iCell = 1:size(tS_session(iMovie).avgTS, 2)
%     
%     curMatTS = squeeze(tS_session(iMovie).matTS_norm(:,iCell,:));
%     avgMatTS = tS_session(iMovie).avgTS_norm(:,iCell);
%     steMatTS = std(curMatTS, [], 2)./sqrt(size(curMatTS, 2)-1);
%     
%     ts(iCell).avgMatTS = avgMatTS;
%     ts(iCell).steMatTS = steMatTS;
%     
%     % figure(100);clf;
%     % subplot(2,1,1)
%     % plot(curMatTS); axis tight
%     % title(sprintf('Cell #%d/%d', iCell, size(tS_session(iMovie).avgTS, 2)))
%     % subplot(2,1,2)
%     % plot(avgMatTS, 'k'); axis tight
%     % line(repmat(1:length(avgMatTS), 2, 1), cat(2, avgMatTS+steMatTS, avgMatTS-steMatTS)', 'Color', 'c')
%     
%     % input('')
% end
% 
% % %% Consistency across trials
% catAvgMatTS = cat(2, ts.avgMatTS); % driven activity, averaged across trials
% catSteMatTS = cat(2, ts.steMatTS);
% 
% [coeff, score, latent, tsquared, explained] = pca(catAvgMatTS');
% 
% k = 4;
% [IDXdfl, C, SUMD] = kmeans(catAvgMatTS', k, 'Distance', 'correlation');
% [sortedIDXdfl, indCelldfl] = sort(IDXdfl);
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
% end
% 
% % Plotting
% fig_summary_DFL = figure;
% set(gcf,  'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1085 750])
% clear sp
% 
% % 1. averaged amplitude for each condition
% figure(fig_summary_DFL);
% sp(1) = subplot('Position', [0.2 0.65 0.75 0.3]);
% imagesc(catAvgMatTS(:, indCelldfl)');
% colormap(hot);
% set(sp(1), 'CLim', [-1 1].*4)
% set(sp(1), 'YTick', find(diff(sortedIDXdfl)>0), 'YTickLabel', [])
% set(sp(1), 'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
% set(sp(1), 'TickDir', 'out')
% box off
% colorbar;
% title(sprintf('%s %s: averaged response to movie %s', nameSubj, dateSession, tS_session(iMovie).idStim))
% ylabel('Cells (sorted)')
% xlabel('Time (s)')
% 
% % 2. Clustering results on 2-d PC space
% figure(fig_summary_DFL);
% sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% for iType = 1:k
%         plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
%         hold on;
% end
% legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Location', 'best')
% xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
% ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
% set(gca, 'TickDir', 'out')
% box off
% axis square
% title('Clustering based on movie response on PC space')
% 
% % 3. Clustering results on imaging field of view
% figure(fig_summary_DFL);
% sp(3) = subplot('Position', [0.55 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% imagesc(imgFOV); 
% colormap(sp(3), gray);
% hold on;
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
%     end
% end
% axis off
% 
% % % for each run
% % fig_map = figure;
% % for iRun = 1:length(tSeries_DFL)
% % matTSnorm = zscore(tSeries_DFL(iRun).C_raw');
% % 
% % k = 5;
% % [IDX, C, SUMD] = kmeans(matTSnorm', k, 'Distance', 'correlation');
% % [sortedIDX, indCell] = sort(IDX);
% % 
% % 
% % clear indCell_sort
% % for iType = 1:k
% %     indCell_sort{iType} = indCell(sortedIDX==iType);
% % end
% % 
% % cMap_sort = hsv(k);
% % figure(fig_map);
% % subplot(1,length(tSeries_DFL),iRun);
% % imagesc(imgFOV); colormap(gray);
% % axis off
% % for iType = 1:k
% %     for iC = 1:size(indCell_sort{iType}, 1)
% %         figure(fig_map);
% %         hold on;
% %         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
% %     end
% % end
% % title(sprintf('DFL Run #%d', iRun))
% % 
% % end
% 
% % end
% 
% % if flagSavePPTX
% %     % save figures
% %     addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
% %     addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
% %     
% %     fname_pptx = sprintf('%s_ClusteringBPMDFL', nameSubj); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
% %     exportFigsToPPTX(fname_pptx);
% %     
% % %     switch lower(nameSubj)
% % %         case 'tabla'
% % %             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
% % %         case 'max'
% % %             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
% % %     end
% %     movefile('./*.pptx', dirFig);
% % end
% 
% % end
% 
% 
% 


