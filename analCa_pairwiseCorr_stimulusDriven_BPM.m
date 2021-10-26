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

flagSavePPTX = 1;

nameSubj ='Tabla'; %'Max'; %

% for iSubj = 1:length(setNameSubj)
% setNameSubj{iSubj}; 

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
%     for iSession = 1:nSession
dateSession = setDateSession{iSession}; %  '20191113'; %
% dateSession =setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_session_stim')

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

%% BPM: averaged response for each image
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_session_stim')
% load(fullfile(dirSave, 'pairwiseCorr.mat'), 'resultsCov');

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


%% DFL: stimulus-evoked responses
iMovie = 1;
matTS = tS_session(iMovie).avgTS_norm;

[matR_DFL] = corr(matTS, 'Type', 'spearman');
[rr, cc, vectR_DFL] = find(triu(matR_DFL, 1));
critAmpR_high = prctile(vectR_DFL, 99.5); %prctile(vectR, 0.5); %prctile(vectR, 99.5); %0.8; %prctile(vectR, 99.5);
critAmpR_low = prctile(vectR_DFL, 0.5); %prctile(vectR, 0.5); %

locHighR = find(vectR_DFL > critAmpR_high);
setPairHighR = cat(2, rr(locHighR), cc(locHighR));
locLowR = find(vectR_DFL < critAmpR_low);
setPairLowR = cat(2, rr(locLowR), cc(locLowR));

imgFOV = neuron.Cn.*neuron.PNR;
cMapType = [0 1 0; 1 0 1; 0 1 1]; %assign green to RS, magenta to BPM, cyan to DFL

iType = 3; % DFL
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

% plot the amplitude profile
fig_selectivity = figure;
set(fig_selectivity, 'Color', 'w');
for iPair = 1:size(setPairHighR, 1)
    figure(fig_selectivity); cla;
    plot(matTS(:, setPairHighR(iPair,:)))
    legend(num2str(setPairHighR(iPair,1)), num2str(setPairHighR(iPair,2)))
    axis tight
    title(sprintf('%s Session %d/: DFL: Pair #%d [Cell %d %d]: r = %2.2f', nameSubj, iSession, iPair, setPairHighR(iPair,1), setPairHighR(iPair,2), vectR_DFL(locHighR(iPair))))
    input('');
end


%% session-by-session data save benchmarking
resultsCorr = struct([]);
        critHigh = 0.005; %0.99;
        for iType = 1:length(matTS_sm) % for diffrent conditions (RS, BPM, DFL)
            switch iType
                case 1
                    nameCond = 'Resting State';
                case 2
                    nameCond = 'Flashing images';
                case 3
                    nameCond = 'Continuous videos';
            end
            
            resultsCorr(iType).nameCond = nameCond;
            
            for iRun = 1:length(matTS_sm(iType).ts)
                clear matR vectR sortedVectR indCellSorted critHighR
                [matR] = corr(matTS_sm(iType).ts{iRun}', 'Type', 'spearman');
                [rr, cc, vectR] = find(triu(matR, 1));
                [sortedVectR, indCellSorted] = sort(vectR);
                critHighR = sortedVectR(round(length(vectR)*(1-critHigh)));
                
                resultsCorr(iType).matR{iRun} = matR;
                resultsCorr(iType).vectR{iRun} = vectR;
                resultsCorr(iType).sortedVectR{iRun} = sortedVectR;
                resultsCorr(iType).indCellSorted{iRun} = indCellSorted;
                resultsCorr(iType).meanR{iRun} = mean(vectR);
                resultsCorr(iType).medianR{iRun} = median(vectR);
                resultsCorr(iType).critHigh = critHigh;
                resultsCorr(iType).critHighR{iRun} = critHighR;
                
            end
        end
        
        resultsCov(iSession).nameSubj = nameSubj;
        resultsCov(iSession).dateSession = dateSession;
        resultsCov(iSession).win_sm = win_sm;
        resultsCov(iSession).matTS_sm = matTS_sm;
        resultsCov(iSession).resultsCorr = resultsCorr;
        
        % save the data
        save(fullfile(dirSave, 'pairwiseCorr.mat'), 'resultsCov');


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


%%%%%%%%%%%%%%%%%%% bis hier
%%
% General distribution of pair-wise corelation
fig_corrDistribution = figure;
set(fig_corrDistribution, 'Color', 'w', 'Position', [100 100 1400 375]);

for iType = 1:length(resultsCorr)
    
    iRun = 1; % look at 1st run here
    
    figure(fig_corrDistribution);
    subplot(1, length(resultsCorr), iType);
    yyaxis left
    hist(resultsCorr(iType).sortedVectR{iRun}, 50);
    line([resultsCorr(iType).medianR{iRun} resultsCorr(iType).medianR{iRun}], ylim)
    text(resultsCorr(iType).medianR{iRun}, max(ylim)-10, ...
        sprintf('median r = %2.2f', resultsCorr(iType).medianR{iRun}), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
    yyaxis right
    plot(resultsCorr(iType).sortedVectR{iRun}, 1:length(resultsCorr(iType).sortedVectR{iRun}));
    line([resultsCorr(iType).critHighR{iRun} resultsCorr(iType).critHighR{iRun}], ylim)
    text(resultsCorr(iType).critHighR{iRun}, max(ylim), ...
        {sprintf('r = %2.2f', resultsCorr(iType).critHighR{iRun}), sprintf('(at highest %2.2f%%)', resultsCorr(iType).critHigh*100)}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    xlabel('Pairwise correlation (rho)')
    title(resultsCorr(iType).nameCond)
    
    if iType == 2
        title({sprintf('%s %s', nameSubj, dateSession), resultsCorr(iType).nameCond})
    end
    
end

% Correlation for each condition
setPairCond = nchoosek(1:3, 2);
iRun = 1; % look at 1st run here
locHighRS = find(resultsCorr(1).sortedVectR{iRun} > resultsCorr(1).critHighR{iRun});
indCell_highRS = resultsCorr(1).indCellSorted{iRun}(locHighRS);
cMap_corrRS = jet(length(indCell_highRS));

fig_pairCorrCond = figure;
set(fig_pairCorrCond, 'Color', 'w', 'Position', [300 300 1140 340]);

for iPair = 1:length(setPairCond)
    
    iRun = 1;
    
    figure(fig_pairCorrCond);
    subplot(1, length(resultsCorr), iPair);
    plot(resultsCorr(setPairCond(iPair, 1)).vectR{iRun}, resultsCorr(setPairCond(iPair, 2)).vectR{iRun}, 'k.'); hold on;
    scatter(resultsCorr(setPairCond(iPair, 1)).vectR{iRun}(indCell_highRS), resultsCorr(setPairCond(iPair, 2)).vectR{iRun}(indCell_highRS), ...
        20, cMap_corrRS, 'filled');
    text(-0.9, 0.9, sprintf('r = %2.2f', corr(resultsCorr(setPairCond(iPair, 1)).vectR{iRun}, resultsCorr(setPairCond(iPair, 2)).vectR{iRun}, ...
        'type', 'spearman')), 'Color', 'k')
    text(-0.9, 0.8, sprintf('r = %2.2f', corr(resultsCorr(setPairCond(iPair, 1)).vectR{iRun}(indCell_highRS), resultsCorr(setPairCond(iPair, 2)).vectR{iRun}(indCell_highRS), ...
        'type', 'spearman')), 'Color', 'm')
    xlabel(sprintf('Corr. in %s', resultsCorr(setPairCond(iPair, 1)).nameCond));
    ylabel(sprintf('Corr. in %s', resultsCorr(setPairCond(iPair, 2)).nameCond));
    
end
h = findobj(gcf, 'type', 'axes');
set(h(:), 'Box', 'off', 'TickDir', 'out', 'XLim', [-1 1], 'YLim', [-1 1]);
title(h(2), sprintf('%s %s: Pairwise correlation in each condition (corr in RS > %2.2f in color)', nameSubj, dateSession,  resultsCorr(1).critHighR{iRun}))


% Pairs of neurons show high correlation under each condition
imgFOV = neuron.Cn.*neuron.PNR;
cMapType = [0 1 0; 1 0 1; 0 1 1]; %assign green to RS, magenta to BPM, cyan to DFL

fig_highCorrPairs = figure;
sizeFig{1} = [1640 410];
sizeFig{2} = [1450 485];
set(fig_highCorrPairs, 'Color', 'w', 'Position', [100 100 sizeFig{iSubj}]);

for iType = 1:length(resultsCorr)
    
    iRun = 1; % look at 1st run here
    
    locHighR = find(resultsCorr(iType).vectR{iRun} > resultsCorr(iType).critHighR{iRun});
    setPairHighR = cat(2, rr(locHighR), cc(locHighR));
    cMap_highR = repmat(linspace(1,0.7, length(locHighR))', 1, 3).*cMapType(iType,:);
    
    figure(fig_highCorrPairs);
    subplot(1, length(resultsCorr), iType); %
    imagesc(imgFOV); colormap(gray);
    hold on
    plot(center(:,2), center(:,1), 'w.'); hold on;
    % set(gca, 'YDir', 'reverse'); hold on;
    for iPair = 1:size(setPairHighR, 1)
        plot(center(setPairHighR(iPair, :), 2), center(setPairHighR(iPair, :), 1), 'o-', ...
            'Color', cMap_highR(iPair,:), 'MarkerFaceColor', cMap_highR(iPair,:), 'MarkerEdgeColor', 'none')
        hold on;
    end
    xlabel(resultsCorr(iType).nameCond)
end
h = findobj(fig_highCorrPairs, 'type', 'axes');
set(h, 'XTick', [], 'YTick', []);
title(h(2), sprintf('%s %s: pairs of neurons showing high correlation in each condition (top %2.2f%%)', nameSubj, dateSession, resultsCorr(1).critHigh*100))



%%




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
