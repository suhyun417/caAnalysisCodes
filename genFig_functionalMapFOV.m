% genFig_functionalMapFOV.m
%
% Map out different properties of cells onto the imaging field of view
% 2020/02/05 SHP

clear all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

setNameSubj = {'Tabla', 'Max'};
flagSavePPTX = 0; %1;


% % for iSubj = 1:2
% 
% nameSubj = setNameSubj{iSubj}; %'Tabla';
% 
% % get session info
% [infoSession, opts] = readInfoSession(nameSubj);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);

close all;

nameSubj = 'Tabla';
dateSession = '20191113';

% for iSession = 1:nSession
% dateSession = setDateSession{iSession}; %'20191113'; %setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data
addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));

% get the contours and image field of view
% neuron_b = neuron.batches{1}.neuron;
thr = 0.3; % the lower the smaller (more centralized) the contour
Coor = neuron.get_contours(thr);
imgFOV = neuron.Cn.*neuron.PNR;

% % draw all the contours
% neuron_b.show_contours([], [], imgFOV, 'true');


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
figure(fig_summary);
sp(1) = subplot('Position', [0.2 0.65 0.75 0.3]);
imagesc(matAmpCellCond(indCell,:)');
colormap(jet);
set(sp(1), 'CLim', [-1 1].*1.5)
if iSubj == 2
    set(sp(1), 'CLim', [-1 1].*1)
end
set(sp(1), 'XTick', find(diff(sortedIDX)>0), 'XTickLabel', [])
set(sp(1), 'YTick', 1:length(setCond), 'YTickLabel', setCondName)
set(sp(1), 'TickDir', 'out')
box off
colorbar;
title(sprintf('%s %s: average response amplitude for each category', nameSubj, dateSession))
xlabel('Cells')
ylabel('Image category')

% 2. Clustering results on 2-d PC space
figure(fig_summary);
sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
cMap_sort = hsv(k);

for iType = 1:k
        plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
        hold on;
end
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Location', 'best')
xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
set(gca, 'TickDir', 'out')
box off
axis square
title('Clustering based on category selectivity on PC space')

% 3. Clustering results on imaging field of view
figure(fig_summary);
sp(3) = subplot('Position', [0.55 0.1 0.4 0.4]);
cMap_sort = hsv(k);

imagesc(imgFOV); 
colormap(sp(3), gray);
hold on;
for iType = 1:k
    for iC = 1:size(indCell_sort{iType}, 1)
        plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
    end
end
axis off


%% From movie: covariation
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))

ts = struct([]);
iMovie = 1;
for iCell = 1:size(tS_session(iMovie).avgTS, 2)
    
    curMatTS = squeeze(tS_session(iMovie).matTS_norm(:,iCell,:));
    avgMatTS = tS_session(iMovie).avgTS_norm(:,iCell);
    steMatTS = std(curMatTS, [], 2)./sqrt(size(curMatTS, 2)-1);
    
    ts(iCell).avgMatTS = avgMatTS;
    ts(iCell).steMatTS = steMatTS;
    
    % figure(100);clf;
    % subplot(2,1,1)
    % plot(curMatTS); axis tight
    % title(sprintf('Cell #%d/%d', iCell, size(tS_session(iMovie).avgTS, 2)))
    % subplot(2,1,2)
    % plot(avgMatTS, 'k'); axis tight
    % line(repmat(1:length(avgMatTS), 2, 1), cat(2, avgMatTS+steMatTS, avgMatTS-steMatTS)', 'Color', 'c')
    
    % input('')
end

% %% Consistency across trials : Well screw this consistency stuff. 
catAvgMatTS = cat(2, ts.avgMatTS); % driven activity, averaged across trials
catSteMatTS = cat(2, ts.steMatTS);

[coeff, score, latent, tsquared, explained] = pca(catAvgMatTS');

k = 4;
[IDXdfl, C, SUMD] = kmeans(catAvgMatTS', k, 'Distance', 'correlation');
[sortedIDXdfl, indCelldfl] = sort(IDXdfl);
clear indCell_sort
for iType = 1:k
    indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
end

% Plotting
fig_summary_DFL = figure;
set(gcf,  'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1085 750])
clear sp

% 1. averaged amplitude for each condition
figure(fig_summary_DFL);
sp(1) = subplot('Position', [0.2 0.65 0.75 0.3]);
imagesc(catAvgMatTS(:, indCelldfl)');
colormap(hot);
set(sp(1), 'CLim', [-1 1].*4)
set(sp(1), 'YTick', find(diff(sortedIDXdfl)>0), 'YTickLabel', [])
set(sp(1), 'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
set(sp(1), 'TickDir', 'out')
box off
colorbar;
title(sprintf('%s %s: averaged response to movie %s', nameSubj, dateSession, tS_session(iMovie).idStim))
ylabel('Cells (sorted)')
xlabel('Time (s)')

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

imagesc(imgFOV); 
colormap(sp(3), gray);
hold on;
for iType = 1:k
    for iC = 1:size(indCell_sort{iType}, 1)
        plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
    end
end
axis off

% % for each run
% fig_map = figure;
% for iRun = 1:length(tSeries_DFL)
% matTSnorm = zscore(tSeries_DFL(iRun).C_raw');
% 
% k = 5;
% [IDX, C, SUMD] = kmeans(matTSnorm', k, 'Distance', 'correlation');
% [sortedIDX, indCell] = sort(IDX);
% 
% 
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCell(sortedIDX==iType);
% end
% 
% cMap_sort = hsv(k);
% figure(fig_map);
% subplot(1,length(tSeries_DFL),iRun);
% imagesc(imgFOV); colormap(gray);
% axis off
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         figure(fig_map);
%         hold on;
%         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
%     end
% end
% title(sprintf('DFL Run #%d', iRun))
% 
% end

% end

if flagSavePPTX
    % save figures
    addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
    addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
    
    fname_pptx = sprintf('%s_ClusteringBPMDFL', nameSubj); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
    exportFigsToPPTX(fname_pptx);
    
%     switch lower(nameSubj)
%         case 'tabla'
%             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
%         case 'max'
%             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
%     end
    movefile('./*.pptx', dirFig);
end

% end

% figure
% [center] = neuron_b.estCenter();
% imagesc(imgFOV)
% hold on
% plot(center(:,2), center(:, 1), 'r.')
% text(center(:,2)+1, center(:,1), num2str([1:95]'), 'Color', 'w');

