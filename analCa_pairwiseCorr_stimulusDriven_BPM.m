% analCa_pairwiseCorr_stimulusDriven_BPM.m
%
% Script to compute correlations based on their stimulus selectivity
% averaged response across trials from flashing image (BPM) runs 
% 2020/03/11 SHP

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

    
nameSubj ='Tabla'; % setNameSubj{iSubj}; %'Max'; %

switch lower(nameSubj)
    case 'tabla'
        dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
    case 'max'
        dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
end

% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[cc, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = cc(2:end); % 1st one is always empty
nSession = length(setDateSession);

close all;

iSession = 2; %1;
dateSession = setDateSession{iSession};
% dateSession = '20191113'; %setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data and compute center coordinates of cells
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

[center] = neuron.estCenter();
center = center(validIndCell, :);

%% pairwise correlation based on stimulus-evoked response
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_session_stim')
load(fullfile(dirSave, 'pairwiseCorr.mat'), 'resultsCov');

matAmpCellStim_a = reshape(cat(1, tS_session_stim.avgAmp_norm), size(tS_session_stim));
matAmpCellStim_b = reshape(cat(1, tS_session_stim.avgAmp_b_norm), size(tS_session_stim));

matAmpCellStim = (matAmpCellStim_a - matAmpCellStim_b)';

% pair-wise correlation 
[matR] = corr(matAmpCellStim, 'Type', 'spearman');
[rr, cc, vectR] = find(triu(matR, 1));
critAmpR_high = prctile(vectR, 99.5); %prctile(vectR, 0.5); %prctile(vectR, 99.5); %0.8; %prctile(vectR, 99.5);
critAmpR_low = prctile(vectR, 0.5); %prctile(vectR, 0.5); %

locHighR = find(vectR > critAmpR_high);
setPairHighR = cat(2, rr(locHighR), cc(locHighR));
locLowR = find(vectR < critAmpR_low);
setPairLowR = cat(2, rr(locLowR), cc(locLowR));


% % plot the amplitude profile
% fig_selectivity = figure;
% set(fig_selectivity, 'Color', 'w');
% for iPair = 1:size(setPairHighR, 1)
%     figure(fig_selectivity); cla;
%     plot(matAmpCellStim(:, setPairHighR(iPair,:)))
%     legend(num2str(setPairHighR(iPair,1)), num2str(setPairHighR(iPair,2)))
%     input('');
% end


imgFOV = neuron.Cn.*neuron.PNR; 
cMapType = [0 1 0; 1 0 1; 0 1 1]; %assign green to RS, magenta to BPM, cyan to DFL

iType = 2; % BPM
cMap_highR = repmat(linspace(1,0.7, length(locHighR))', 1, 3).*cMapType(iType,:);        

fig_highCorrPairs = figure;
% sizeFig{1} = [1640 410];
% sizeFig{2} = [1450 485];
set(fig_highCorrPairs, 'Color', 'w') %, 'Position', [100 100 sizeFig{iSubj}]);

figure(fig_highCorrPairs);
% subplot(1, length(resultsCorr), iType); %
imagesc(imgFOV); colormap(gray);
hold on
plot(center(:,2), center(:,1), 'w.'); hold on;
% set(gca, 'YDir', 'reverse'); hold on;
for iPair = 1:size(setPairHighR, 1)
    plot(center(setPairHighR(iPair, :), 2), center(setPairHighR(iPair, :), 1), 'o-', ...
        'Color', cMap_highR(iPair,:), 'MarkerFaceColor', cMap_highR(iPair,:), 'MarkerEdgeColor', 'none')
    hold on;
end


% cf. single trial across neuron correlation
resultsCorr = resultsCov(iSession).resultsCorr;

iType = 3;
iRun = 1; % look at 1st run here
locHighR = find(resultsCorr(iType).vectR{iRun} > resultsCorr(iType).critHighR{iRun});
[rr, cc, vectR] = find(triu(resultsCorr(iType).matR{iRun}, 1));
setPairHighR = cat(2, rr(locHighR), cc(locHighR));
fig_corrTS = figure;
for iPair = 1:size(setPairHighR,1)
figure(fig_corrTS); clf;
plot(resultsCov(iSession).matTS_sm(iType).ts{iRun}(setPairHighR(iPair,:), :)');
axis tight
title(sprintf('%s Session %d/: %s: Pair #%d [Cell %d %d]: r = %2.2f', nameSubj, iSession, resultsCorr(iType).nameCond, iPair, setPairHighR(iPair,1), setPairHighR(iPair,2), resultsCorr(iType).vectR{iRun}(locHighR(iPair))))
input('')
end
