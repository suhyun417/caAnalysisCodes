% analCa_covariation.m

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
        
        %% Pairwise correlation
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts.mat', nameSubj, dateSession), 'tSeries_BPM')
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/RS_ts.mat', nameSubj, dateSession), 'tSeries_RS')
        tS_session_DFL = tS_session;
        clear tS_session
        
        % compute smoothed time series
        matTS_sm = struct([]); % 1 for resting state, 2 for BPM, 3 for DFL
        win_sm = 20; % 2 sec window smoothing
        for iRun = 1:length(tSeries_RS)
            matTS_sm(1).ts{iRun} = smoothdata(zscore(tSeries_RS(iRun).C_raw(validIndCell, :), [], 2), 2, 'gaussian', win_sm);
        end
        for iRun = 1:length(tSeries_BPM)
            matTS_sm(2).ts{iRun} = smoothdata(zscore(tSeries_BPM(iRun).C_raw(validIndCell, :), [], 2), 2, 'gaussian', win_sm);
        end
        countRunDFL = 0;
        for iMovie = 1:2
            for iRun = 1:size(tS_session_DFL(iMovie).matTS_norm, 3)
                countRunDFL = countRunDFL+1;
                matTS_sm(3).ts{countRunDFL} = smoothdata(tS_session_DFL(iMovie).matTS_norm(:,validIndCell, iRun), 'gaussian', win_sm)';
            end
        end
        for iType = 1:length(matTS_sm)
            switch iType
                case 1
                    nameCond = 'Resting State';
                case 2
                    nameCond = 'Flashing images';
                case 3
                    nameCond = 'Continuous videos';
            end
            
            matTS_sm(iType).nameCond = nameCond;
        end
        % matTS_sm.tsDFL = smoothdata(tS_session_DFL(1).avgTS_norm(:,validIndCell), 'gaussian', win_sm)';
        
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
        
        % Pairs that show high correlation in resting state but low in visual experiments
        % iType = 2;
        critLow = 0.05; %
        
        iRun = 1; % look at 1st run here
        locHighRS = find(resultsCorr(1).sortedVectR{iRun} > resultsCorr(1).critHighR{iRun});
        indCell_highRS = resultsCorr(1).indCellSorted{iRun}(locHighRS);
        
        fig_highCorrRS_lowCorrVis = figure;
        set(fig_highCorrRS_lowCorrVis, 'Color', 'w', 'Position', [300 300 1075 395]);
        if iSubj==2
            set(fig_highCorrRS_lowCorrVis, 'Color', 'w', 'Position', [300 300 950 486]);
        end
        
        for iType = 2:3 % only for visual experiments
            
            locLowVisR = find(abs(resultsCorr(iType).vectR{iRun}(indCell_highRS))<critLow);
            setPair_lowVisR = cat(2, rr(indCell_highRS(locLowVisR)), cc(indCell_highRS(locLowVisR)));
            cMap_lowR = repmat(linspace(1, 0.9, length(locLowVisR))', 1, 3).*cMapType(iType,:);
            
            figure(fig_highCorrRS_lowCorrVis);
            subplot(1,2,iType-1)
            imagesc(imgFOV); colormap(gray);
            hold on
            plot(center(:,2), center(:,1), 'w.'); hold on;
            % set(gca, 'YDir', 'reverse'); hold on;
            for iPair = 1:size(setPair_lowVisR, 1)
                plot(center(setPair_lowVisR(iPair, :), 2), center(setPair_lowVisR(iPair, :), 1), 'o-', ...
                    'Color', cMap_lowR(iPair,:), 'MarkerFaceColor', cMap_lowR(iPair,:), 'MarkerEdgeColor', 'none')
                hold on;
            end
            set(gca, 'XTick', [], 'YTick', [])
            xlabel(sprintf('Pairs with |r|<%2.2f during %s', critLow, resultsCorr(iType).nameCond))
        end
        h = findobj(gcf, 'type', 'axes');
        % set(h(:), 'Box', 'off', 'TickDir', 'out', 'XLim', [-1 1], 'YLim', [-1 1]);
        % title(h(2), sprintf('%s %s: pairs of neurons with r > %2.2f during RS', nameSubj, dateSession, resultsCorr(1).critHighR))
        mtit(sprintf('%s %s: pairs of neurons show high r (> %2.2f) during resting but show low r in visual experiments', nameSubj, dateSession, resultsCorr(1).critHighR{iRun}))
        
    end
    
    
    if flagSavePPTX
        % save figures
        addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
        addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
        
        fname_pptx = sprintf('figs_pairwiseCorr_%s', nameSubj); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
        exportFigsToPPTX_SHP(fname_pptx);
        
        movefile(sprintf('./*%s.pptx', nameSub), dirSave);
    end
    
end





% lowDFL = [10 11 14 18 19 21 22 29 37 55];
% figure
% cat(2, r(ind(lowDFL)), c(ind(lowDFL)))
%
% plot(center([28 33], 2), center([28 33], 1), 'o', 'MarkerFaceColor', cMap(1,:)); hold on
% plot(center([31    34], 2), center([31    34], 1), 'o', 'MarkerFaceColor', cMap(2,:))
% plot(center([8    45], 2), center([8    45], 1), 'o', 'MarkerFaceColor', cMap(3,:))
% plot(center([60    61], 2), center([60    61], 1), 'o', 'MarkerFaceColor', cMap(4,:))
% plot(center([51    63], 2), center([51    63], 1), 'o', 'MarkerFaceColor', cMap(5,:))
% plot(center([55    70], 2), center([55    70], 1), 'o', 'MarkerFaceColor', cMap(6,:))
% plot(center([52 58 66], 2), center([52 58 66], 1), 'o', 'MarkerFaceColor', cMap(7, :))
% plot(center([72 81 109], 2), center([72 81 109], 1), 'o', 'MarkerFaceColor', cMap(8, :))
% xlim([1 size(neuron.Cn, 2)])
% ylim([1 size(neuron.Cn, 1)])
% set(gca, 'YDir', 'reverse');
% lowFlash = [3 6 42 46 47 49 50 57];
% cat(2, r(ind(lowFlash)), c(ind(lowFlash)))
% plot(center([12 16 17], 2), center([12 16 17], 1), '^', 'MarkerFaceColor', cMap(end,:))
% plot(center([83 84 87 95], 2), center([83 84 87 95], 1), '^', 'MarkerFaceColor', cMap(end-1,:))
% plot(center([26 102], 2), center([26 102], 1), '^', 'MarkerFaceColor', cMap(end-2,:))
% plot(center([22 103], 2), center([22 103], 1), '^', 'MarkerFaceColor', cMap(end-3,:))
% plot(center([81   110], 2), center([81   110], 1), '^', 'MarkerFaceColor', cMap(end-4,:))


% %% Resting state
% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/RS_ts.mat', nameSubj, dateSession))
%
% matTS_norm = zscore(tSeries_RS(1).C_raw(validIndCell, :), [], 2);
% win_sm = 20; % 2 sec window smoothing
% matTS_norm_smooth = smoothdata(matTS_norm, 2, 'gaussian', win_sm);
%
% matTS = matTS_norm_smooth; %matTS_norm; %matTS_norm_smooth;
%
% [coeff, score, latent, tsquared, explained] = pca(matTS);
% [sorted_score1, indCell_score1] = sort(score(:,1), 'descend');
% figure
% imagesc(matTS(indCell_score1, :))
% colormap(hot)
% set(gca, 'CLim', [0 1].*5)
% figure
% imagesc(neuron.PNR);
% colormap(gray);
% hold on
% crit_score =20;
% plot(center(indCell_score1(abs(sorted_score1)<crit_score), 2), center(indCell_score1(abs(sorted_score1)<crit_score), 1), 'w.', 'MarkerSize', 12);
% plot(center(indCell_score1(sorted_score1<-crit_score), 2), center(indCell_score1(sorted_score1<-crit_score), 1), 'c.', 'MarkerSize', 12);
% plot(center(indCell_score1(sorted_score1>crit_score), 2), center(indCell_score1(sorted_score1>crit_score), 1), 'm.', 'MarkerSize', 12);
%
% [matR] = corr(matTS', 'Type', 'spearman');
% [rr, cc, vectR] = find(triu(matR, 1));
% figure; hist(vectR, 50);
% ind = find(vectR>0.6);
% ind2d = sub2ind(size(matR), rr(ind), cc(ind));
%
% % check the spatial location of correlated neurons
% cMap = cool(length(ind));
% [sortedR, sortedInd] = sort(matR(ind2d));
%
% figCenter = figure;
% for iPair = 1:length(ind)
%     plot(center(rr(ind(sortedInd(iPair))), 2), center(rr(ind(sortedInd(iPair))), 1), 'o', 'MarkerFaceColor', cMap(iPair, :), 'MarkerEdgeColor', cMap(iPair, :))
%     hold on;
%     plot(center(cc(ind(sortedInd(iPair))), 2), center(cc(ind(sortedInd(iPair))), 1), 'o', 'MarkerFaceColor', cMap(iPair, :), 'MarkerEdgeColor', cMap(iPair, :))
% end
% set(gca, 'YDir', 'reverse');
% xlim([1 size(neuron.Cn, 2)])
% ylim([1 size(neuron.Cn, 1)])
% setCell = unique(cat(1, rr(ind), cc(ind)));
% text(center(setCell, 2)+1, center(setCell, 1), num2str(setCell), 'Color', 'k')
%
%
%
% % % quick clustering
% % k = 4;
% % [IDX, C, SUMD] = kmeans(matTS, k);
% % [sortedIDX, indCell] = sort(IDX);
% % figure; imagesc(matTS(indCell, :)); % quick check
% % colormap(hot)
% % set(gca, 'CLim', [0 5])
% %
% % % check the spatial clustering
% % for iType = 1:k
% %     indCell_sort{iType} = indCell(sortedIDX==iType);
% % end
% %
% % imgFOV = neuron.Cn;
% %
% % cMap_sort = hsv(k);
% % fig_map = figure;
% % imagesc(imgFOV); colormap(gray);
% % axis off
% % for iType = 1:k
% %     for iC = 1:size(indCell_sort{iType}, 1)
% %         figure(fig_map);
% %         hold on;
% %         plot(center(indCell_sort{iType}(iC, 1), 2), center(indCell_sort{iType}(iC, 1), 1), '.', 'Color', cMap_sort(iType, :), 'MarkerSize', 10);
% %     end
% % end
%
%
% %% driven correlation across neurons during movie and PC
% for iSession = 1:nSession
%     dateSession = setDateSession{iSession};
%
%     dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%     dirPreproc = fullfile(dirProcdata_session, '_preproc');
%
%     load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
%     tS_session_DFL = tS_session;
%     clear tS_session
%
%     ts_DFL = struct([]);
%     for iMovie = 1:2
%         matTS = tS_session_DFL(iMovie).avgTS_norm;
%         win_sm = 20; % 2 sec window smoothing
%         matTS_sm = smoothdata(matTS, 'gaussian', win_sm);
%         ts_DFL(iMovie).matTS = matTS;
%         ts_DFL(iMovie).matTS_sm = matTS_sm;
%     end
%
%     % [matR_DFL] = corr(tS_DFL(iMovie).matTS_sm, 'Type', 'spearman');
%     % [r_dfl, c_dfl, vectR_DFL] = find(triu(matR_DFL, 1));
%     % figure; hist(vectR_DFL, 50);
%     % ind = find(vectR_DFL>0.6);
%     % ind2d = sub2ind(size(matR), r_dfl(ind), c_dfl(ind));
%
%     fig_PC = figure;
%     set(fig_PC, 'Color', 'w', 'Position', [94   596   854   520]);
%     for iMovie = 1:2
%         [coeff, score, latent, tsquared, explained] = pca(ts_DFL(iMovie).matTS_sm');
%
%         subplot(2,1,iMovie)
%         plot(coeff(:,1:3))
%         axis tight
%         legend(sprintf('PC1 (%d%%)', round(explained(1))), sprintf('PC2 (%d%%)', round(explained(2))), sprintf('PC3 (%d%%)', round(explained(3))))
%         title(sprintf('%s %s: Eigenvectors of movie driven activity for movie %d (2s smooth)', nameSubj, dateSession, iMovie))
%         set(gca, 'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
%         xlabel('Time (s)')
%     end
%     print(fig_PC, sprintf(fullfile(dirFig, 'movieDriven2secSmooth_PCA_PCs_%s_%s'), nameSubj, dateSession), '-depsc')
% end
%
% %% Dynamic face localizer (movie)
% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
%
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
% % %% Consistency across trials : Well screw this consistency stuff.
% catAvgMatTS = cat(2, ts.avgMatTS);
% catSteMatTS = cat(2, ts.steMatTS);
% %
% % % clear hist*
% % % bin = -2:0.5:7;
% % % for iCell = 1:length(ts)
% % % [N,EDGES] = histcounts(ts(iCell).avgMatTS, bin);
% % % [Nste,EDGESste] = histcounts(ts(iCell).steMatTS,bin);
% % % histAvg(:,iCell) = N';
% % % histSte(:,iCell) = Nste';
% % % end
% %
% % clear indCell_consistency
% % avgSteMatTS = mean(catSteMatTS);
% % [sortSte, indCell] = sort(avgSteMatTS, 'descend');
% % indCell_consistency{1} = indCell(1:round(length(indCell)*0.1))';
% % indCell_consistency{2} = indCell(end-round(length(indCell)*0.1)+1:end)';
% %
% % neuron_b = neuron.batches{1}.neuron;
% % thr = 0.3; % the lower the smaller (more centralized) the contour
% % Coor = neuron_b.get_contours(thr);
% % imgFOV = neuron_b.Cn.*neuron_b.PNR;
% %
% % % neuron_b.show_contours([], neuron.ids(indCell(1:round(length(indCell)*0.1))), neuron_b.PNR.*neuron_b.Cn, 'true');
% %
% % cMap_sort = cool(2); % blue to pink
% % fig_map = figure;
% % imagesc(imgFOV); colormap(gray);
% % axis off
% % for iType = 1:2
% %     for iC = 1:size(indCell_consistency{iType}, 1)
% %         figure(fig_map);
% %         hold on;
% %         plot(Coor{indCell_consistency{iType}(iC, 1)}(1,:), Coor{indCell_consistency{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
% %     end
% % end
% % title(sprintf('%s %s', nameSubj, dateSession))
% %
% % addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
% % addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
% %
% % fname_pptx = 'DFL_acrossTrialConsistency'; % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
% % exportFigsToPPTX(fname_pptx);
% %
% % movefile(sprintf('./%s.pptx', fname_pptx), dirFig);
%
%
% %% visually driven activity (i.e. averaged across trial)
% catAvgMatTS = cat(2, ts.avgMatTS);
%
% k = 5;
% [IDX, C, SUMD] = kmeans(catAvgMatTS', k, 'Distance', 'correlation');
% [sortedIDX, indCell] = sort(IDX);
% figure;
%
%
% [matR] = corr(catAvgMatTS, 'Type', 'spearman');
%
% %% Load the cell time series
% load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
% load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
%
% iRun = 1;
% matTSWhole_BPM = tSeries_BPM(iRun).C_raw';
% matTSWhole_BPM_norm = zscore(matTSWhole_BPM);
% [matR] = corr(matTSWhole_BPM_norm, 'Type', 'spearman');
%
% matTSWhole_DFL = tSeries_DFL(iRun).C_raw';
% matTSWhole_DFL_norm = zscore(matTSWhole_DFL);
% [matR] = corr(matTSWhole_DFL_norm, 'Type', 'spearman');
%
%
% % quick clustering
% k = 5;
% [IDX, C, SUMD] = kmeans(matAmpCellStim, k);
% [sortedIDX, indCell] = sort(IDX);
% figure; imagesc(matAmpCellStim(indCell, :)'); % quick check
% colormap(jet)
% set(gca, 'CLim', [-3 3])
%
% % check the spatial clustering
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
%
%
% %% quick clustering based on each run's time series
% fig_map = figure;
% for iRun = 1:length(tSeries_BPM)
% matTSnorm = zscore(tSeries_BPM(iRun).C_raw');
%
% k = 5;
% [IDX, C, SUMD] = kmeans(matTSnorm', k);
% [sortedIDX, indCell] = sort(IDX);
%
%
%
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCell(sortedIDX==iType);
% end
%
% cMap_sort = hsv(k);
% figure(fig_map);
% subplot(1,k,iRun);
% imagesc(imgFOV); colormap(gray);
% axis off
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         figure(fig_map);
%         hold on;
%         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
%     end
% end
% title(sprintf('BPM Run #%d', iRun))
%
% end
%
% figure; imagesc(matTSnorm(:,indCell)'); % quick check
% set(gca, 'YTick', find(diff(sortedIDX)>0))
% colormap(jet)
% set(gca, 'CLim', [-1 1].*5)
%
% %% quick clustering based on each run's time series
% fig_map = figure;
% for iRun = 1:length(tSeries_DFL)
% matTSnorm = zscore(tSeries_DFL(iRun).C_raw');
%
% k = 5;
% [IDX, C, SUMD] = kmeans(matTSnorm', k, 'Distance', 'correlation');
% [sortedIDX, indCell] = sort(IDX);
%
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