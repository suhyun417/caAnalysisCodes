% genFig_LNtalk.m

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


nameSubj = 'Tabla';
dateSession = '20191125'; % '20191113';

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
neuron_b = neuron.batches{1}.neuron;
thr = 0.3; % the lower the smaller (more centralized) the contour
Coor = neuron_b.get_contours(thr);
imgFOV = neuron_b.Cn.*neuron_b.PNR;

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