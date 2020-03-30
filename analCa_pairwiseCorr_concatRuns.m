% analCa_pairwiseCorr_concatRuns.m
%
% Compute the pairwise correlation of concatenated time series for each BPM
% and DFL conditions
% modified from "analCa_pairwiseCorr.m"
% 2020/03/23 SHP

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

flagSavePPTX = 0; %1;

setNameSubj = {'Tabla', 'Max'};

for iSubj = 1:length(setNameSubj)
    
    nameSubj = setNameSubj{iSubj}; %'Max'; %'Tabla';
    
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
    
    for iSession = 1:nSession
        % iSession = 1;
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
        
        %% Concatenate the timeseries for each condition
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_session')
        tS_session_BPM_trial = tS_session.tS_trial;
        clear tS_session;
        %         load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts.mat', nameSubj, dateSession), 'tSeries_BPM')
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
        tS_session_DFL = tS_session;
        clear tS_session
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/RS_ts.mat', nameSubj, dateSession), 'tSeries_RS')
        
        % compute smoothed time series
        matTS_sm = struct([]); % 1 for resting state, 2 for BPM, 3 for DFL
        win_sm = 20; % 2 sec window smoothing
        
        % Resting state
        tTS_RS = [];
        tTS_RS = cat(2, tSeries_RS.C_raw);
        matTS_sm(1).ts = smoothdata(zscore(tTS_RS(validIndCell, :), [], 2), 2, 'gaussian', win_sm);
        
        % BPM
        tTS_BPM = [];
        for iCell = 1:size(tS_session_BPM_trial, 1)
            tTS_BPM(iCell, :) = cat(1, tS_session_BPM_trial(iCell, :).matTS_norm)';
        end
        matTS_sm(2).ts = smoothdata(tTS_BPM(validIndCell, :), 2, 'gaussian', win_sm);        
%         tTS_BPM = [];
%         for iCell = 1:size(tS_session_BPM_trial, 1)
%             tTS_BPM(iCell, :) = cat(1, tS_session_BPM_trial(iCell, :).matTS)';
%         end
%         matTS_sm(2).ts = smoothdata(zscore(tTS_BPM(validIndCell, :), [], 2), 2, 'gaussian', win_sm);
        
        % DFL
        tTS_DFL = [];
        tTS_DFL = reshape(permute(cat(3, tS_session_DFL.matTS_norm), [2 1 3]), size(tS_session_DFL(1).matTS_norm, 2), []);
        matTS_sm(3).ts = smoothdata(tTS_DFL(validIndCell, :), 2, 'gaussian', win_sm);
%         tTS_DFL = [];
%         tTS_DFL = reshape(permute(cat(3, tS_session_DFL.matTS), [2 1 3]), size(tS_session_DFL(1).matTS, 2), []);
%         matTS_sm(3).ts = smoothdata(zscore(tTS_DFL(validIndCell, :), [], 2), 2, 'gaussian', win_sm);
        
        clear tTS*
        
        
        %% Compute pairwise correlation
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
            
            clear matR vectR sortedVectR indCellSorted critHighR
            [matR] = corr(matTS_sm(iType).ts', 'Type', 'spearman');
            [rr, cc, vectR] = find(triu(matR, 1));
            [sortedVectR, indCellSorted] = sort(vectR);
            critHighR = sortedVectR(round(length(vectR)*(1-critHigh)));
            
            resultsCorr(iType).matR = matR;
            resultsCorr(iType).vectR = vectR;
            resultsCorr(iType).sortedVectR = sortedVectR;
            resultsCorr(iType).indCellSorted = indCellSorted;
            resultsCorr(iType).meanR = mean(vectR);
            resultsCorr(iType).medianR = median(vectR);
            resultsCorr(iType).critHigh = critHigh;
            resultsCorr(iType).critHighR = critHighR;
            
        end
        
        resultsCov(iSession).nameSubj = nameSubj;
        resultsCov(iSession).dateSession = dateSession;
        resultsCov(iSession).win_sm = win_sm;
        resultsCov(iSession).matTS_sm = matTS_sm;
        resultsCov(iSession).center = center;
        resultsCov(iSession).neuron.Cn = neuron.Cn;
        resultsCov(iSession).neuron.PNR = neuron.PNR;
        resultsCov(iSession).resultsCorr = resultsCorr;
        
        % save the data
        save(fullfile(dirSave, 'pairwiseCorr_20sSM_concat.mat'), 'resultsCov');
        
% %         % General distribution of pair-wise corelation
% %         fig_corrDistribution = figure;
% %         set(fig_corrDistribution, 'Color', 'w', 'Position', [100 100 1400 375]);
% %         
% %         for iType = 1:length(resultsCorr)
% %             
% %             figure(fig_corrDistribution);
% %             subplot(1, length(resultsCorr), iType);
% %             yyaxis left
% %             hist(resultsCorr(iType).sortedVectR, 50);
% %             line([resultsCorr(iType).medianR resultsCorr(iType).medianR], ylim)
% %             text(resultsCorr(iType).medianR, max(ylim)-10, ...
% %                 sprintf('median r = %2.2f', resultsCorr(iType).medianR), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
% %             yyaxis right
% %             plot(resultsCorr(iType).sortedVectR, 1:length(resultsCorr(iType).sortedVectR));
% %             line([resultsCorr(iType).critHighR resultsCorr(iType).critHighR], ylim)
% %             text(resultsCorr(iType).critHighR, max(ylim), ...
% %                 {sprintf('r = %2.2f', resultsCorr(iType).critHighR), sprintf('(at highest %2.2f%%)', resultsCorr(iType).critHigh*100)}, ...
% %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
% %             xlabel('Pairwise correlation (rho)')
% %             title(resultsCorr(iType).nameCond)
% %             
% %             if iType == 2
% %                 title({sprintf('%s %s', nameSubj, dateSession), resultsCorr(iType).nameCond})
% %             end
% %             
% %         end
% %         h = findobj(gcf, 'type', 'axes');
% %         set(h(:), 'Box', 'off', 'TickDir', 'out', 'XLim', [-1 1]); %, 'YLim', [-1 1]);
% %         
% %         % Correlation for each condition
% %         setPairCond = nchoosek(1:3, 2);
% %         iRun = 1; % look at 1st run here
% %         locHighRS = find(resultsCorr(1).sortedVectR > resultsCorr(1).critHighR);
% %         indCell_highRS = resultsCorr(1).indCellSorted(locHighRS);
% %         cMap_corrRS = jet(length(indCell_highRS));
% %         
% %         fig_pairCorrCond = figure;
% %         set(fig_pairCorrCond, 'Color', 'w', 'Position', [300 300 1140 340]);
% %         
% %         for iPair = 1:length(setPairCond)
% %             
% %             figure(fig_pairCorrCond);
% %             subplot(1, length(resultsCorr), iPair);
% %             plot(resultsCorr(setPairCond(iPair, 1)).vectR, resultsCorr(setPairCond(iPair, 2)).vectR, 'k.'); hold on;
% %             scatter(resultsCorr(setPairCond(iPair, 1)).vectR(indCell_highRS), resultsCorr(setPairCond(iPair, 2)).vectR(indCell_highRS), ...
% %                 20, cMap_corrRS, 'filled');
% %             text(-0.9, 0.9, sprintf('r = %2.2f', corr(resultsCorr(setPairCond(iPair, 1)).vectR, resultsCorr(setPairCond(iPair, 2)).vectR, ...
% %                 'type', 'spearman')), 'Color', 'k')
% %             text(-0.9, 0.8, sprintf('r = %2.2f', corr(resultsCorr(setPairCond(iPair, 1)).vectR(indCell_highRS), resultsCorr(setPairCond(iPair, 2)).vectR(indCell_highRS), ...
% %                 'type', 'spearman')), 'Color', 'm')
% %             xlabel(sprintf('Corr. in %s', resultsCorr(setPairCond(iPair, 1)).nameCond));
% %             ylabel(sprintf('Corr. in %s', resultsCorr(setPairCond(iPair, 2)).nameCond));
% %             
% %         end
% %         h = findobj(gcf, 'type', 'axes');
% %         set(h(:), 'Box', 'off', 'TickDir', 'out', 'XLim', [-1 1], 'YLim', [-1 1]);
% %         title(h(2), sprintf('%s %s: Pairwise correlation in each condition (corr in RS > %2.2f in color)', nameSubj, dateSession,  resultsCorr(1).critHighR))
% %         
% %         
% %         % Pairs of neurons show high correlation under each condition
% %         imgFOV = neuron.Cn.*neuron.PNR;
% %         cMapType = [0 1 0; 1 0 1; 0 1 1]; %assign green to RS, magenta to BPM, cyan to DFL
% %         
% %         fig_highCorrPairs = figure;
% %         sizeFig{1} = [1640 410];
% %         sizeFig{2} = [1450 485];
% %         set(fig_highCorrPairs, 'Color', 'w', 'Position', [100 100 sizeFig{iSubj}]);
% %         
% %         for iType = 1:length(resultsCorr)
% %             
% %             locHighR = find(resultsCorr(iType).vectR > resultsCorr(iType).critHighR);
% %             setPairHighR = cat(2, rr(locHighR), cc(locHighR));
% %             cMap_highR = repmat(linspace(1,0.7, length(locHighR))', 1, 3).*cMapType(iType,:);
% %             
% %             figure(fig_highCorrPairs);
% %             subplot(1, length(resultsCorr), iType); %
% %             imagesc(imgFOV); colormap(gray);
% %             hold on
% %             plot(center(:,2), center(:,1), 'w.'); hold on;
% %             % set(gca, 'YDir', 'reverse'); hold on;
% %             for iPair = 1:size(setPairHighR, 1)
% %                 plot(center(setPairHighR(iPair, :), 2), center(setPairHighR(iPair, :), 1), 'o-', ...
% %                     'Color', cMap_highR(iPair,:), 'MarkerFaceColor', cMap_highR(iPair,:), 'MarkerEdgeColor', 'none')
% %                 hold on;
% %             end
% %             xlabel(resultsCorr(iType).nameCond)
% %         end
% %         h = findobj(fig_highCorrPairs, 'type', 'axes');
% %         set(h, 'XTick', [], 'YTick', []);
% %         title(h(2), sprintf('%s %s: pairs of neurons showing high correlation in each condition (top %2.2f%%)', nameSubj, dateSession, resultsCorr(1).critHigh*100))
        
%         % Pairs that show high correlation in resting state but low in visual experiments
%         % iType = 2;
%         critLow = 0.05; %
%         
%         locHighRS = find(resultsCorr(1).sortedVectR > resultsCorr(1).critHighR);
%         indCell_highRS = resultsCorr(1).indCellSorted(locHighRS);
%         
%         fig_highCorrRS_lowCorrVis = figure;
%         set(fig_highCorrRS_lowCorrVis, 'Color', 'w', 'Position', [300 300 1075 395]);
%         if iSubj==2
%             set(fig_highCorrRS_lowCorrVis, 'Color', 'w', 'Position', [300 300 950 486]);
%         end
%         
%         for iType = 2:3 % only for visual experiments
%             
%             locLowVisR = find(abs(resultsCorr(iType).vectR(indCell_highRS))<critLow);
%             setPair_lowVisR = cat(2, rr(indCell_highRS(locLowVisR)), cc(indCell_highRS(locLowVisR)));
%             cMap_lowR = repmat(linspace(1, 0.9, length(locLowVisR))', 1, 3).*cMapType(iType,:);
%             
%             figure(fig_highCorrRS_lowCorrVis);
%             subplot(1,2,iType-1)
%             imagesc(imgFOV); colormap(gray);
%             hold on
%             plot(center(:,2), center(:,1), 'w.'); hold on;
%             % set(gca, 'YDir', 'reverse'); hold on;
%             for iPair = 1:size(setPair_lowVisR, 1)
%                 plot(center(setPair_lowVisR(iPair, :), 2), center(setPair_lowVisR(iPair, :), 1), 'o-', ...
%                     'Color', cMap_lowR(iPair,:), 'MarkerFaceColor', cMap_lowR(iPair,:), 'MarkerEdgeColor', 'none')
%                 hold on;
%             end
%             set(gca, 'XTick', [], 'YTick', [])
%             xlabel(sprintf('Pairs with |r|<%2.2f during %s', critLow, resultsCorr(iType).nameCond))
%         end
%         h = findobj(gcf, 'type', 'axes');
%         % set(h(:), 'Box', 'off', 'TickDir', 'out', 'XLim', [-1 1], 'YLim', [-1 1]);
%         % title(h(2), sprintf('%s %s: pairs of neurons with r > %2.2f during RS', nameSubj, dateSession, resultsCorr(1).critHighR))
%         mtit(sprintf('%s %s: pairs of neurons show high r (> %2.2f) during resting but show low r in visual experiments', nameSubj, dateSession, resultsCorr(1).critHighR))
        
    end
    
    
    if flagSavePPTX
        % save figures
        addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
        addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
        
        fname_pptx = sprintf('figs_pairwiseCorr_%s_concatenated', nameSubj); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
        exportFigsToPPTX_SHP(fname_pptx);
        
        movefile(sprintf('./*%s*.pptx', nameSubj), dirSave);
    end
    
end

%% Checking the time series
nameSubj = 'Max'; %'Tabla';

switch lower(nameSubj)
    case 'tabla'
        dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
    case 'max'
        dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
end

% Load covariation results
load(fullfile(dirSave, 'pairwiseCorr_20sSM_concat.mat'), 'resultsCov');

iSession = 1; %3;
resultsCorr = resultsCov(iSession).resultsCorr;
center = resultsCov(iSession).center;
imgFOV = neuron.Cn.*neuron.PNR;

iType = 3;
locHighR = find(resultsCorr(iType).vectR > resultsCorr(iType).critHighR);
[rr, cc, vectR] = find(triu(resultsCorr(iType).matR, 1));
setPairHighR = cat(2, rr(locHighR), cc(locHighR));

fig_corrTS = figure;
set(fig_corrTS, 'Color', 'w', 'Position', [60 560 1210 300])
cMapType = [0 1 0; 1 0 1; 0 1 1]; %assign green to RS, magenta to BPM, cyan to DFL
cMap_highR = repmat(linspace(1,0.7, length(locHighR))', 1, 3).*cMapType(iType,:);
for iPair = 1:size(setPairHighR,1)
    figure(fig_corrTS); clf;
    AX(1) = subplot('Position', [0.01 0.1 0.15 0.8]);
    imagesc(imgFOV);
    set(AX(1), 'Colormap', gray)
    imagesc(imgFOV); colormap(gray); axis off;
    hold on
    plot(center(:,2), center(:,1), 'w.'); hold on;
    plot(center(setPairHighR(iPair, :), 2), center(setPairHighR(iPair, :), 1), '.-', ...
                    'Color', cMap_highR(iPair,:), 'MarkerFaceColor', cMap_highR(iPair,:), 'MarkerEdgeColor', 'none');
                text(center(setPairHighR(iPair, :), 2)+1, center(setPairHighR(iPair, :), 1), num2str(setPairHighR(iPair,:)'), 'Color', cMap_highR(iPair,:))
    AX(2) = subplot('Position', [0.22 0.15 0.75 0.75]);
    plot(resultsCov(iSession).matTS_sm(iType).ts(setPairHighR(iPair,:), :)');
    axis tight
    AX(2).XTickLabel = cellstr(num2str(str2double(AX(2).XTickLabel)./10));
    xlabel(AX(2), 'Time (s)')
    ylabel(AX(2), 'Norm. Resp. (sd)')
    title(AX(2), sprintf('%s Session %d: %s: Pair #%d [Cell %d %d]: r = %2.2f', nameSubj, iSession, resultsCorr(iType).nameCond, iPair, setPairHighR(iPair,1), setPairHighR(iPair,2), resultsCorr(iType).vectR(locHighR(iPair))))
    input('')
end
