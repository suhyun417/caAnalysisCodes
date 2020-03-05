% analCa_pairwiseCorr_summary.m

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

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% flagSavePPTX = 1;

setNameSubj = {'Tabla', 'Max'};

for iSubj = 1:length(setNameSubj)
    
    nameSubj = setNameSubj{iSubj}; %'Max'; %'Tabla';
    
    switch lower(nameSubj)
        case 'tabla'
            dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
        case 'max'
            dirSave = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
    end
    
    %% Load covariation results
    load(fullfile(dirSave, 'pairwiseCorr.mat'), 'resultsCov');
    
    pairHighCorr = struct([]);
    for iSession = 1:length(resultsCov)
        resultsCorr = resultsCov(iSession).resultsCorr;
        center = resultsCov(iSession).center;
        
        for iType = 1:length(resultsCorr)
            
            for iRun = 1:length(resultsCorr(iType).matR)
                [rr, cc, vectR] = find(triu(resultsCorr(iType).matR{iRun}, 1));
                
                clear locHighR setPairHighR
                locHighR = find(resultsCorr(iType).vectR{iRun} > resultsCorr(iType).critHighR{iRun});
                setPairHighR = cat(2, rr(locHighR), cc(locHighR));
                d = sqrt((center(setPairHighR(:,1), 1) - center(setPairHighR(:,2), 1)).^2 + (center(setPairHighR(:,1), 2) - center(setPairHighR(:,2), 2)).^2);
                
                pairHighCorr(iSession, iType).setPairHighR{iRun} = setPairHighR;
                pairHighCorr(iSession, iType).setD{iRun} = d;
                pairHighCorr(iSession, iType).meanD{iRun} = mean(d);
                pairHighCorr(iSession, iType).medianD{iRun} = median(d);
                pairHighCorr(iSession, iType).steD{iRun} = std(d)/sqrt(length(d)-1);
            end
            
            pairHighCorr(iSession, iType).grandMeanD = mean(cat(1, pairHighCorr(iSession, iType).meanD{:}));
            pairHighCorr(iSession, iType).grandSteD = std(cat(1, pairHighCorr(iSession, iType).setD{:}))/sqrt(length(cat(1, pairHighCorr(iSession, iType).setD{:}))-1);
            
        end
        
    end
    
    matGrandMeanD = reshape(cat(1, pairHighCorr.grandMeanD), size(pairHighCorr));
    matGrandSteD = reshape(cat(1, pairHighCorr.grandSteD), size(pairHighCorr));
    
    figure;
    set(gcf, 'Color', 'w')
    cMap_session = cool(size(pairHighCorr,1));
    set(gca, 'ColorOrder', cMap_session)
    hold on;
    for iSession = 1:size(pairHighCorr, 1)
        line([1:3;1:3], [matGrandMeanD(iSession,:)-matGrandSteD(iSession,:); matGrandMeanD(iSession,:)+matGrandSteD(iSession,:)], 'Color', cMap_session(iSession, :), 'LineWidth', 2)
        hold on;
    end
    plot(1:3, matGrandMeanD, 'o-', 'MarkerFaceColor', 'w', 'LineWidth', 2, 'MarkerSize', 10)
    xlim([0.5 3.5])
    set(gca, 'Box', 'off', 'TickDir', 'out')
    set(gca, 'XTick', 1:3)
    set(gca, 'XTickLabel', {resultsCorr.nameCond})
    ylabel('Euclidean distance between cells')
    
    title(sprintf('%s', nameSubj))
    print(gcf, sprintf(fullfile(dirFig, 'distance_highCorrNeuronPairs_%s'), nameSubj), '-depsc')
    
end

for iSession = 1:length(resultsCov)
        resultsCorr = resultsCov(iSession).resultsCorr;
        center = resultsCov(iSession).center;

        for iType = 1:3
            clear B I
            [B, I] = sort(resultsCov(iSession).resultsCorr(iType).matR{1}, 'descend');
        
            sortedMatD = []; sortedMatR = [];
            for iCell = 1:size(B, 1)
                clear center_pairs
                curCell = ones(size(B, 1)-1, 1).*iCell;
                %             center_cell1 = center(curCell, :);
                %             center_cell2 = center(I(2:end, iCell), :);
                center_pairs(:,:,1) = center(curCell,:);
                center_pairs(:,:,2) = center(I(2:end, iCell), :);
                d = sum(diff(center_pairs, 1, 3).^2, 2);
                sortedMatD(:,iCell) = d;
            end
            sortedMatR = B(2:end, :);
            
            corrD(iSession, iType).sortedMatD = sortedMatD;
            corrD(iSession, iType).sortedMatR = sortedMatR;
        end
end

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
    
nc = 30;
edges = -1:0.1:1;
cMap = jet(length(edges)-1);
figSingleCell = figure;
for iCell = 1:size(resultsCov(iSession).resultsCorr(iType).matR{1}, 1)
    figure(figSingleCell); clf;
    for iType = 1:3
        
        [B, I] = sort(resultsCov(iSession).resultsCorr(iType).matR{1}, 'ascend');
        indCorr = discretize(B(:,iCell), edges);
        
        subplot(1, 3, iType);        
        set(gca, 'YDir', 'reverse'); hold on;
        scatter(center(I(:, iCell), 2), center(I(:, iCell), 1), 30, cMap(indCorr, :), 'fill'); %jet(nc), 'fill');
        text(center(iCell, 2)+1, center(iCell, 1), num2str(iCell), 'Color', 'k');
        title(sprintf('%s Session %d/%d: Cell %d', nameSubj, iSession, length(resultsCov), iCell))
        
    end
    
    input('')
end




scatter(center(I(:, iCell), 2), center(I(:, iCell), 1), 30, cMap(indCorr, :), 'fill');
    
    
    
    
    