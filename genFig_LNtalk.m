% genFig_LNtalk.m

clear all;

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

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %1; %2; %1;

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

% aligned cells TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')

%%
% nameSubj = 'Tabla';
% dateSession = '20191125'; % '20191113';

cMap_line = cool(length(setDateSession));
% fig_contours = figure;
for iSession = 1:length(setDateSession) %;
dateSession = setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data
addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));

% get the contours and image field of view
% thr = 0.3; % the lower the smaller (more centralized) the contour
% Coor = neuron.get_contours(thr);
% imgFOV = neuron.Cn.*neuron.PNR;
% 
% % draw all the contours
% neuron.show_contours([], [], imgFOV, 'true');

% draw contours for each session
figure;
subplot('Position', [0 0 1 1]);
imagesc(neuron.Cn.*neuron.PNR); colormap(gray);
hold on;

thr = 0.6; %0.2;
cellColor = cMap_line(iSession, :); %'m'; %'c'; 
widthContour = 1;
[d1,d2] = size(neuron.Cn);
indCellValid_session = cellIDAcrossDay(~isnan(cellIDAcrossDay(:,iSession)), iSession);

CC = cell(length(indCellValid_session),1);
CR = cell(length(indCellValid_session),2);
% cmap_cell = cool(size(neuron.A, 2));
for iCC = 1:length(indCellValid_session)
    i = indCellValid_session(iCC); %sortedIndCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
truesize;

print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_cMapCool', ...
    nameSubj, FOV_ID, thr*10, iSession)), '-depsc');

end



% %% Stacked cells across days: playing with examples
% indCellValid = find(cat(1, cellTS.nTrial1_total)>8);
% 
% % tempATrial = cat(2, cellPix(indCellValid).repPix);
% % tempATrial(~isnan(tempATrial)) = 10;
% tempA = cat(2, cellPix.repPix);
% tempA(~isnan(tempA)) = 1;
% tempA(:, indCellValid) = tempA(:, indCellValid).*10;
% 
% imgCells = sum(tempA, 2, 'omitnan');
% imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));
% 
% figure;
% set(gcf, 'Color', 'w')
% imagesc(imgCells_2d)
% colormap(turbo)
% 
% % A_temp = full(reshape(neuron.A(:,i),d1,d2));
% 
% figure;
% set(gcf, 'Color', 'w')
% cmap_cell = colormap(hsv(length(cellPix)));
% for iCell = 1:length(cellPix)
%     plot(cellPix(iCell).contourCell{1}(1,1:end), cellPix(iCell).contourCell{1}(2,1:end), ...
%         'Color', cmap_cell(iCell, :), 'linewidth', 1); hold on;
% %     text(cellPix(iCell).contourCell{1}(1,end), cellPix(iCell).contourCell{1}(2,end), num2str(iCell), ...
% %         'color', 'k')
% end
% set(gca, 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])


%% Visualization of cell contours across days
cellColor = cool(length(setDateSession));

fig_acSession = figure;
set(fig_acSession, 'Color', 'w')

for iSession = 1:length(setDateSession)
    % iSession = 1;
    dateSession = setDateSession{iSession}; %'20191113'; %setDateSession{iSession};
    
    dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    
    %%
    addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
    cnmfe_setup;
    d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
    
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
    
    % validIndCell = [];
    % validIndCell(:,1) = 1:length(neuron.ids);
    % if strcmpi(nameSubj, 'max')
    %     load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
    %     validIndCell = indCell.validCell;
    % end
    
    validIndCell = cellIDAcrossDay(~isnan(cellIDAcrossDay(:,iSession)), iSession);
    
    
    % color cells to indicate sessions
    thr = 0.3;
    widthContour = 1; %1.5;
    Coor = neuron.get_contours(thr);
    infoCells(iSession).coor_0p3 = Coor;
    
    figure(fig_acSession);
    subplot('Position', [0 0 1 1])
    for i = 1:length(validIndCell) %size(Coor, 1)
        %         cont = medfilt1(Coor{i}')';
        cont = Coor{validIndCell(i)};
        if size(cont,2) > 1 % "shifts" values are in the order of image matrix dimensions: first dimension is vertical ("y") and second dimension is horizontal ("x")
            plot(cont(1,1:end)+shifts(iSession, 2), cont(2,1:end)+shifts(iSession, 1), 'Color', cellColor(iSession, :), 'linewidth', widthContour); hold on;
        end
    end
end
if iSubj == 1
    set(gca, 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])
    axis equal
    axis off
    set(fig_acSession, 'Position', [1200 1140 size(infoCells(1).imgFOV,2)+40 size(infoCells(1).imgFOV,1)+40])
    print(fig_acSession, fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_AcSession_cMapCool_padding40', ...
    nameSubj, FOV_ID, thr*10)), '-depsc');
elseif iSubj ==2
    set(gca, 'YDir', 'reverse', 'XLim', [0 size(infoCells(1).imgFOV,2)], 'YLim', [0 size(infoCells(1).imgFOV,1)])
    axis equal
    axis off
    set(fig_acSession, 'Position', [1200 1140 size(infoCells(1).imgFOV,2) size(infoCells(1).imgFOV,1)])
    print(fig_acSession, fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_AcSession_cMapCool', ...
    nameSubj, FOV_ID, thr*10)), '-depsc');
end



tempA = cat(2, cellPix.repPix);
tempA(~isnan(tempA)) = 1;
imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));
aa = cat(2, zeros(size(imgCells_2d,1), 20), imgCells_2d, zeros(size(imgCells_2d,1),20));
aa = cat(1, zeros(20, size(infoCells(1).imgFOV,2)+40), aa, zeros(20, size(infoCells(1).imgFOV,2)+40));
figure;
subplot('Position', [0 0 1 1])
imagesc(aa, [0 1]);
colormap(gray)
axis off
set(gcf, 'Position', [1200 1140 size(infoCells(1).imgFOV,2)+40 size(infoCells(1).imgFOV,1)+40])
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_validCellContourf_AllCellsAcSession_gray_padding40', ...
    nameSubj, FOV_ID)), '-depsc');

tempA = cat(2, cellPix.repPix);
tempA(~isnan(tempA)) = 1;
imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));
figure;
subplot('Position', [0 0 1 1])
imagesc(imgCells_2d, [0 1]);
colormap(gray)
set(gcf, 'Position', [1200 1140 size(infoCells(1).imgFOV,2) size(infoCells(1).imgFOV,1)])
axis off
print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_validCellContourf_AllCellsAcSession_gray', ...
    nameSubj, FOV_ID)), '-depsc');

% for thr = [0.2:0.1:0.5]
%     locbad=[]; nrows=[]; ncols=[];
% Coor = neuron.get_contours(thr); 
% [nrows, ncols] = cellfun(@size, Coor);
% locbad = find(ncols<2); % the ones that get_contours couldn't get a reasonable localized contour at this threshold
% length(locbad)
% end

figure;
for iC = 1:length(locbad)    
    i = locbad(iC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', 'b'); hold on;
    end
    set(gca, 'YDir', 'reverse', 'XLim', [0 d2], 'YLim', [0 d1])
    input('')
end




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

% thr = 0.3; % the lower the smaller (more centralized) the contour
Coor = neuron.get_contours(thr);
imgFOV = neuron.Cn.*neuron.PNR;


figure;
[center] = neuron.estCenter();
imagesc(imgFOV); colormap(gray);
hold on
plot(center(:,2), center(:, 1), 'r.');
text(center(:,2)+1, center(:,1), num2str([1:length(neuron.ids)]'), 'Color', 'w');


% draw contours
figure;
subplot('Position', [0 0 1 1]);
imagesc(neuron.Cn.*neuron.PNR); colormap(gray);
hold on;

thr = 0.6; %0.2;
% cellColor = [1 1 1];
linecolor = 'm'; %'c'; 
widthContour = 1;
[d1,d2] = size(neuron.Cn);
indCellValid_session = cellIDAcrossDay(~isnan(cellIDAcrossDay(:,iSession)), iSession);

CC = cell(length(indCellValid_session),1);
CR = cell(length(indCellValid_session),2);
% cmap_cell = cool(size(neuron.A, 2));
for iCC = 1:length(indCellValid_session)
    i = indCellValid_session(iCC); %sortedIndCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', linecolor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off

title(sprintf('%s: %s', nameSubj, dateSession))



%% LN talk 2020 Feb
%%
condName_BPM = {'human face', 'marmoset face', 	'marmoset body', 'scene', 'non familiar object', 'hands and catcher',...
    'phase scrambled', 'space scrambled', 'grating', 'random dot motion'};

setCond = [1 2 5 6 10]; %cat(2, 11:15, 21:25, 51
nImage = 5;

catCondMat = cat(1, stimTiming_BPM.condMat);
curSetCond = unique(catCondMat(:,1));
if sum(ismember(setCond, curSetCond))~=5
    setCond = curSetCond;
end

setCondName = condName_BPM(setCond);

%%  for each category
matAmpCellStim = reshape(cat(1, tS_session_stim.avgAmp_norm), size(tS_session_stim));
catStimCondSession = floor(cat(1, tS_session_stim(1,:).idStim)./10);
matAmpCellCond = [];
for iCond = 1:length(setCond)
    idCond = setCond(iCond);
    ind = find(catStimCondSession == idCond);
    matAmpCellCond(:,iCond) = mean(matAmpCellStim(:,ind), 2);
end

[coeff, score, latent, tsquared, explained] = pca(matAmpCellCond);
cumsum(explained);

k = 4;
[IDX, C, SUMD] = kmeans(matAmpCellCond, k);
[sortedIDX, indCell] = sort(IDX);
clear indCell_sort
for iType = 1:k
indCell_sort{iType} = indCell(sortedIDX==iType);
end


fig_summary = figure;
set(gcf, 'Position', [100 100 1085 750])
clear sp

% 1. averaged amplitude for each condition
figure;
sp(1) = gca; %subplot('Position', [0.2 0.65 0.75 0.3]);
imagesc(matAmpCellCond(indCell,:)');
colormap(jet);
set(sp(1), 'CLim', [-1 1].*1.5)
% if iSubj == 2
%     set(sp(1), 'CLim', [-1 1].*1)
% end
set(sp(1), 'XTick', find(diff(sortedIDX)>0), 'XTickLabel', [])
set(sp(1), 'YTick', 1:length(setCond), 'YTickLabel', setCondName)
set(sp(1), 'TickDir', 'out')
box off
colorbar;
title(sprintf('%s %s: average response amplitude for each category', nameSubj, dateSession))
xlabel('Cells')
ylabel('Image category')

% % 2. Clustering results on 2-d PC space
% figure(fig_summary);
% sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% for iType = 1:k
%         plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
%         hold on;
% end
% legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Location', 'best')
% xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
% ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
% set(gca, 'TickDir', 'out')
% box off
% axis square
% title('Clustering based on category selectivity on PC space')

% 3. Clustering results on imaging field of view
figure;
cMap_sort = hsv(k);

imagesc(imgFOV); 
colormap(gray);
hold on;
for iType = 1:k
%     iType = setClust(iK);
    for iC = 1:size(indCell_sort{iType}, 1)
        plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
    end
end
axis off

% all the neurons and 
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(imgFOV); hold on;
colormap(gray);
cMap_all = spring(length(idNeuron));
for iCell = 1:length(Coor)
    plot(Coor{idNeuron(iCell)}(1,:), Coor{idNeuron(iCell)}(2,:), '.', 'Color', cMap_all(iCell,:));
    hold on;
end

% %% summary responses for each item
% matAmpCellStim = reshape(cat(1, tS_session_stim.avgAmp), size(tS_session_stim));
% figure;
% set(gcf, 'Color', 'w', 'Position', [680         602        1080         376])
% imagesc(matAmpCellStim')
% colormap(jet)
% set(gca, 'CLim', [-1 1])
% xlabel('Cells')
% ylabel('Stimulus')
% %         condName = {infoTrial.infoStim([1:6:25]).nameCondition};
% set(gca, 'YTickLabel', setCondName)
% 
% % Quick clustering
% k = 5;
% [IDX, C, SUMD] = kmeans(matAmpCellStim, k);
% [sortedIDX, indCell] = sort(IDX);
% figure; imagesc(matAmpCellStim(indCell, :)'); % quick check
% colormap(jet)
% set(gca, 'CLim', [-1 1])
% 
% % check the spatial clustering
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCell(sortedIDX==iType);
% end
% 
% cMap_sort = hsv(k);
% fig_map = figure;
% imagesc(imgFOV); colormap(gray);
% axis off
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         figure(fig_map);
%         hold on;
%         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
%     end
% end