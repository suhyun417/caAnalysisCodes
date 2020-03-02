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

flagSavePPTX = 1;

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
        
        movefile(sprintf('./*%s.pptx', nameSubj), dirSave);
    end
    
end
