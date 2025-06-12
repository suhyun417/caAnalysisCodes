% genFig_responseConsistency.m
%
% 2024/05/13 SHP
% visualize the consistency across trials/sessions
% focus on finding the most effective way to show consistency of
%   - example neurons, populations
%   - over different timescale
%   - over different aspects, e.g. selectivity, fluctuations etc.
% likely the figure 4 of manuscript

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

%%
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

% aligned cells movie TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')



%% trial-to-trial heatmap of movie 1 responses for each neuron
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 

setCell_ID = [2, 22, 141, 14, 31, 34, 181, 75]; % Indices of example neurons (out of 337 cells)

for iCell = 1:length(setCell_ID)
    idCell = setCell_ID(iCell);
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 100 850 170])
    subplot('Position', [0 0 1 1])
    imagesc(zscore(cellTS(idCell).matTS_movie1')');
    colormap("hot")
    set(gca, 'Box', 'off', 'CLim', [-1 9])
    axis off
    print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_respHeatmap_alltrials_cellID%03d_337_cLim-19', nameSubj, FOV_ID, idCell)), '-depsc')
    print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_respHeatmap_alltrials_cellID%03d_337_cLim-19', nameSubj, FOV_ID, idCell)), '-r300', '-dtiff')
end


%% Spatial map
% get an average center coords
centerCoords = [];
for iC = 1:length(cellPix)
    curMean = mean(cellPix(iC).centerCell, 1);
    centerCoords(iC, 1) = curMean(1);
    centerCoords(iC, 2) = curMean(2);
end


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

setCell_ID = [2, 22, 141, 14, 31, 34, 181, 75]; % Indices of example neurons (out of 337 cells)

figure;
set(gcf, 'Color', 'w')
% cMap_sort = hsv(k);
% cMap_sort(2,:) = [206 182 49]./255
% cMap_sort = [215 25 28; 253 174 97; 255 255 191; 171 217 233; 44 123 182]./255;
% cMap_sort = [208 28 139; 241 182 218; 247 247 247; 184 225 134; 77 172 38]./255;
[d1 d2] = size(infoCells(1).imgFOV);
% imagesc(zeros(size(infoCells(1).imgFOV)));
% colormap(gray)
% imagesc(imgFOV);
% colormap(sp(3), gray);
% hold on;
% for iType = 1:k
for iC = 1:size(cellPix, 2)
Coor = cellPix(iC).contourCell{1};
p = plot(Coor(1,:), Coor(2,:), '.', 'Color', ones(1,3).*0.75); hold on;
if ismember(iC, setCell_ID)
%     locC = find(setCell_ID==iC);
    p.Color = [0 0 0];
end
% text(Coor(1,end), Coor(2,end), num2str(indCellValid(indCell_sort{iType}(iC, 1))), ...
%         'color', 'k')
end
% end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20]) %, 'Color', 'k', 'Box', 'off')
axis equal
axis off
line([300 300-(100/3.125)], [290 290], 'color', 'k', 'LineWidth', 5) % scale bar for 100 um (1px = 3.125um)
set(gca, 'XTick', [], 'YTick', [], 'Box', 'off')
setCellID_str = sprintf('%d_', setCell_ID);
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_allCellgray0p75_setCellID%s_black',nameSubj, FOV_ID, setCellID_str)), '-depsc')
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_allCellgray0p75_setCellID%s_black',nameSubj, FOV_ID, setCellID_str)), '-r300', '-dtiff')


%% session TS examples
figtemp = figure;
set(gcf, 'color', 'w', 'Position', [401    35   378   862])

% tempSetC = find(matDaysAcRegistration>25);
% for iC = 1:length(tempSetC)

setCell_ID = [2, 22, 141, 14, 31, 34, 181, 75]; % Indices of example neurons (out of 337 cells)

for iCell = 1:length(setCell_ID)
idCell = setCell_ID(iCell); %104; %315; %tempSetC(iC); %34; %23;
curMatTS_movie1 = [];
curMatTS_movie1 = zscore(cellTS(idCell).matTS_movie1')';


setT = []; avgT=[];
setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);
for iSS = 1:length(setT)
    avgT(:,iSS) = mean(curMatTS_movie1(setT(iSS,1):setT(iSS,2),:),1)'; %mean(cellTS(idCell).matTS_movie1(setT(iSS,1):setT(iSS,2),:),1)';
end


stringNameSession= cat(1, setDateSession{cellTS(idCell).idAcrossSession(:,1)}); % get date
cMap = cool(length(setT));

figtemp = figure;
set(gcf, 'color', 'w', 'Position', [401    35   378   862])
figure(figtemp); cla;
for iSS = 1:length(setT)
figure(figtemp);
ydata = avgT(:,iSS)+5*(iSS-1);
plot(ydata, 'color', cMap(iSS,:)); %, 'LineWidth', 1.5);
ytick(iSS) = mean(ydata);
hold on;
end
axis tight
set(gca, 'TickDir', 'out', 'box', 'off')
set(gca, 'YTick', ytick, 'YTicklabel', stringNameSession)
set(gca,'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
title(sprintf('%s: Cell ID %d', nameSubj, idCell))
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession_z',nameSubj, FOV_ID, idCell)), '-depsc')
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession_z',nameSubj, FOV_ID, idCell)), '-r300', '-dtiff')


figtemp2 = figure;
set(gcf, 'color', 'w', 'Position', [401    35   556   289])
alpha = 0.2;
figure(figtemp2); cla;
for iSS = 1:length(setT)
figure(figtemp2);
% ydata = avgT(:,iSS); %+50*(iSS-1);
p = plot(avgT(:,iSS)); %, 'color', cMap(iSS,:)); %, 'LineWidth', 1.5);
p.Color = [cMap(iSS,:) alpha];
p.LineWidth = 2;
% ytick(iSS) = mean(ydata);
hold on;
end
axis tight
set(gca, 'TickDir', 'out', 'box', 'off')
% set(gca, 'YTick', ytick, 'YTicklabel', stringNameSession)
set(gca,'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
title(sprintf('%s: Cell ID %d', nameSubj, idCell))
print(figtemp2, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession_overlayZ0p2',nameSubj, FOV_ID, idCell)), '-depsc')
print(figtemp2, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession_overlayZ0p2',nameSubj, FOV_ID, idCell)), '-r300', '-dtiff')

% input('')
end


%% Case of Max
%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %2; %1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%%
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

% aligned cells movie TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')



%% trial-to-trial heatmap of movie 1 responses for each neuron
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % 

setCell_ID = [2, 12, 13, 36, 49, 54, 58, 60]; % Indices of example neurons (out of 171 cells)

for iCell = 1:length(setCell_ID)
    idCell = setCell_ID(iCell);
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 100 850 170])
    subplot('Position', [0 0 1 1])
    imagesc(zscore(cellTS(idCell).matTS_movie1')');
    colormap("hot")
    ttt = get(gca, 'CLim');
%     matLim(iCell,:) = ttt;
    set(gca, 'Box', 'off', 'CLim', [-1 ttt(2)])
    axis off
    print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_respHeatmap_alltrials_cellID%03d_171_cLim-1max', nameSubj, FOV_ID, idCell)), '-depsc')
    print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_respHeatmap_alltrials_cellID%03d_171_cLim-1max', nameSubj, FOV_ID, idCell)), '-r300', '-dtiff')
end


%% Spatial map
% get an average center coords
centerCoords = [];
for iC = 1:length(cellPix)
    curMean = mean(cellPix(iC).centerCell, 1);
    centerCoords(iC, 1) = curMean(1);
    centerCoords(iC, 2) = curMean(2);
end


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

setCell_ID = [2, 12, 13, 36, 49, 54, 58, 60]; % Indices of example neurons (out of 171 cells)

figure;
set(gcf, 'Color', 'w')
% cMap_sort = hsv(k);
% cMap_sort(2,:) = [206 182 49]./255
% cMap_sort = [215 25 28; 253 174 97; 255 255 191; 171 217 233; 44 123 182]./255;
% cMap_sort = [208 28 139; 241 182 218; 247 247 247; 184 225 134; 77 172 38]./255;
[d1 d2] = size(infoCells(1).imgFOV);
% imagesc(zeros(size(infoCells(1).imgFOV)));
% colormap(gray)
% imagesc(imgFOV);
% colormap(sp(3), gray);
% hold on;
% for iType = 1:k
for iC = 1:size(cellPix, 2)
Coor = cellPix(iC).contourCell{1};
p = plot(Coor(1,:), Coor(2,:), '.', 'Color', ones(1,3).*0.75); hold on;
if ismember(iC, setCell_ID)
%     locC = find(setCell_ID==iC);
    p.Color = [0 0 0];
end
% text(Coor(1,end), Coor(2,end), num2str(indCellValid(indCell_sort{iType}(iC, 1))), ...
%         'color', 'k')
end
% end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 d2+20], 'YLim', [0-20 d1+20]) %, 'Color', 'k', 'Box', 'off')
axis equal
axis off
line([0 0], [250 250-(100/3.125)], 'color', 'k', 'LineWidth', 5) % scale bar for 100 um (1px = 3.125um)
set(gca, 'XTick', [], 'YTick', [], 'Box', 'off')
setCellID_str = sprintf('%d_', setCell_ID);
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_allCellgray0p75_setCellID%s_black',nameSubj, FOV_ID, setCellID_str)), '-depsc')
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_allCellgray0p75_setCellID%s_black',nameSubj, FOV_ID, setCellID_str)), '-r300', '-dtiff')






%%%%%%
%% Corelation within and across sessions (obsolete as of 2025/05, as we decided to exclude)
%%%%%%
% flagCell = ~isnan(cellIDAcrossDay);
% 
% setDateSession_datenum = datenum(setDateSession, 'yyyymmdd');
% nDaysBetweenSessions = diff(setDateSession_datenum);
% 
% matDaysAcRegistration = [];
% nSessionRegistered = sum(flagCell, 2);
% for iCell = 1:length(nSessionRegistered)
%     if nSessionRegistered(iCell) > 1
%         locSession = find(flagCell(iCell, :)>0);
%         nDays = setDateSession_datenum(locSession(end)) - setDateSession_datenum(locSession(1));
%         matDaysAcRegistration(iCell,1) = nDays;
%     else
%         matDaysAcRegistration(iCell,1) = 0;
%     end
% end
% 
% %% within session across trial correlation & between session trial-average correlation
% resultsCorr = struct([]);
% for iCell = 1:length(matDaysAcRegistration)
%     idCell = iCell;
%     
%     resultsCorr(iCell).maxDays = 0; % default
%      
%     if nSessionRegistered(iCell) > 1  
% 
%         resultsCorr(iCell).idCell = idCell;
% 
%         locSession = find(flagCell(idCell, :)>0);
%         maxDays = setDateSession_datenum(locSession(end)) - setDateSession_datenum(locSession(1));
%         setDApart = diff(setDateSession_datenum(locSession));
%         resultsCorr(iCell).maxDays = maxDays;
%         resultsCorr(iCell).setDApart = setDApart;
% 
%         setT = [];
%         setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);
%         numT = diff(setT, 1, 2)+1;
% 
%         ssTS = []; bssetrr_setD = [];
%         for iSSS = 1:length(numT)
% 
%             rr = corr(cellTS(idCell).matTS_movie1(setT(iSSS,1):setT(iSSS,2),:)', 'type', 'Spearman');
%             setrr = tril(rr,-1);
%             setrr = setrr(abs(setrr)>0);
%             meanrr = mean(setrr);
% 
%             resultsCorr(iCell).setR{iSSS} = setrr;
%             resultsCorr(iCell).meanR(iSSS,1) = meanrr;
% 
%             % across-trial average TS
%             ssTS(:, iSSS) = mean(cellTS(idCell).matTS_movie1(setT(iSSS,1):setT(iSSS,2),:),1)';   
% 
%             if iSSS < length(numT)
%             bssetrr_setD = cat(1, bssetrr_setD, cumsum(setDApart(iSSS:end)));
%             end
% 
% %             figure(100);
% %             set(gcf, 'Color', 'w', 'Position', [652   823   906   271]);
% %             subplot(2,1,1)
% %             imagesc(cellTS(idCell).matTS_movie1(setT(iSSS,1):setT(iSSS,2),:))
% %             title(sprintf('Cell %d/%d:, Session %d/%d', iCell, length(matDaysAcRegistration), iSSS, length(numT)))
% %             subplot(2,1,2)
% %             plot(cellTS(idCell).matTS_movie1(setT(iSSS,1):setT(iSSS,2),:)')
% %             axis tight
% %             legend
% %             title(sprintf('mean rho = %2.3f', meanrr))
% %             input('')
% 
%         end
%         
%         resultsCorr(iCell).grandMeanR = mean(resultsCorr(iCell).meanR); % within-session across-trial correlation
% 
%         bsrr = corr(ssTS, 'type', 'Spearman');
%         bssetrr = tril(bsrr,-1);
%         bssetrr = bssetrr(abs(bssetrr)>0);
%         bsmeanrr = mean(bssetrr);
% 
%         
% 
%         resultsCorr(iCell).btsn_setR = bssetrr;
%         resultsCorr(iCell).btsn_meanR = bsmeanrr;
%         resultsCorr(iCell).btsn_setR_setD = bssetrr_setD;
% 
%     else
%         continue;
%     end
% 
% end
% 
% % within-session across trial correlations
% setMaxD = cat(1, resultsCorr.maxDays);
% setCellID = cat(1, resultsCorr(setMaxD>0).idCell); % 199 cells 
% catGrandM_ws = cat(1, resultsCorr.grandMeanR);
% catSetR_ws = []; catSetR_ws_cellID = [];
% for iCell = 1:length(setCellID)
%     idCell = setCellID(iCell);
%     catSetR_ws = cat(1, catSetR_ws, resultsCorr(idCell).setR{:}); % all within-session trial pairs (n = 5697)
%     catSetR_ws_cellID = cat(1, catSetR_ws_cellID, repmat(idCell, length(cat(1, resultsCorr(idCell).setR{:})), 1));
% end
% 
% % between session correlations
% catSetR_bs = cat(1, resultsCorr.btsn_setR); % n=2385 (pair-wise sessions)
% catMeanR_bs = cat(1, resultsCorr.btsn_meanR);
% catSetR_bs_days = cat(1, resultsCorr.btsn_setR_setD);
% catSetR_bs_cellID = [];
% for iCell = 1:length(setCellID)
%     idCell = setCellID(iCell);
%     catSetR_bs_cellID = cat(1, catSetR_bs_cellID, repmat(idCell, length(resultsCorr(idCell).btsn_setR), 1));
% end
% 
% % matSetR_ws = cat(2, catSetR_ws, catSetR_ws_cellID);
% % matSetR_bs = cat(2, catSetR_bs, catSetR_bs_cellID, catSetR_bs_days);
%  
% % compare within-session corr and between-session corr for each cell
% figure;
% set(gcf, 'Color', 'w')
% plot(catGrandM_ws, catMeanR_bs, 'o'); hold on;
% axis square
% xlabel('within-session across trial correlation')
% ylabel('across-session correlation')
% 
% % save(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_movieRespCorr.mat',...
% %     nameSubj, FOV_ID, nameSubj, FOV_ID)), 'resultsCorr');
% 
% 
% %% correlation between 1st and last session that are at least 2 weeks apart
% locWeek = find(matDaysAcRegistration>=14);
% length(locWeek)
% 
% % gather averaged TS of 1st session and last session that is at least 2
% % weeks apart
% matAvgTS_1 = []; matAvgTS_2 = []; tempInfoCell = [];
% for iCell = 1:length(locWeek)
%     idCell = locWeek(iCell);
%     locSession = find(flagCell(idCell, :)>0);
%     nDays = setDateSession_datenum(locSession(end)) - setDateSession_datenum(locSession(1));
% 
%     setT = []; 
%     setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);
% 
%     matAvgTS_1(:,iCell) = mean(cellTS(idCell).matTS_movie1(setT(1,1):setT(1,2),:),1)';
%     matAvgTS_2(:,iCell) = mean(cellTS(idCell).matTS_movie1(setT(end,1):setT(end,2),:),1)';
%     tempInfoCell(iCell, 1) = idCell;
%     tempInfoCell(iCell, 2) = nDays;
% 
% end
% 
% [matR] = corr(matAvgTS_1, matAvgTS_2, 'type', 'Spearman');
% setR = diag(matR);
% figure
% hR = histogram(setR);
% figure
% plot(tempInfoCell(:,2), setR, 'o')
% 
% % quick test
% for iCell = 1:length(locWeek)
%     figure(3); cla;
%     plot(matAvgTS_1(:,iCell), 'ro-');
%     hold on
%     plot(matAvgTS_2(:,iCell), 'bo-');
%     input('')
% end
% 
% 
% 
% %% How the across-session correlation is related to overall response level?
% % iSession = 1; %:length(setDateSession) %;
% % dateSession = setDateSession{iSession};
% % 
% % dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% % dirPreproc = fullfile(dirProcdata_session, '_preproc');
% % 
% % % dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
% % 
% % load(sprintf('%s/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', directory.dirProcdata, nameSubj, dateSession))
% % % load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))
% 
% 
% 
% %% In case I need the histogram values
% % figure;
% % h = histogram(sum(flagCell, 2));
% % 
% % plot(cumsum(fliplr(h.Values)), 'o-')
% % hold on
% % text((1:length(h.Values))'-0.25, (cumsum(fliplr(h.Values))+10)', num2str(cumsum(fliplr(h.Values))'))
% % set(gca, 'XTick', 1:length(h.Values), 'XTickLabel', length(h.Values):-1:1)
% % xlabel('Number of detected sessions')
% % ylabel('Cumulative number of cells')
% % set(gca, 'TickDir', 'out', 'Box', 'off')
% % set(gca, 'FontSize', 15, 'LineWidth', 2)
% % % print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_cellRegistration_cumNumOverSessions', nameSubj, FOV_ID)), '-r300', '-depsc')
% 


