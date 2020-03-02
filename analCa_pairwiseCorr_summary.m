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

    
    
    

    

    
    
    
    
    