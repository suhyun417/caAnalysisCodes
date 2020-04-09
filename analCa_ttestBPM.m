% analCa_ttestBPM.m
%
% script to evaluate visually responsive neurons based on responses
% from flashing image (BPM) runs using t-test with Bonferroni correction
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


setNameSubj = {'Tabla', 'Max'};

dirSave{1} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Tabla/FOV1');
dirSave{2} = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/Max/FOV3');

flagPlot = 1;

for iSubj = 1:length(setNameSubj)
    nameSubj = setNameSubj{iSubj}; %'Max'; %'Tabla'; %
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj);
    
    [cc, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = cc(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    resultsTTest = struct([]);
    
    for iSession = 1:length(setDateSession)
        
        dateSession = setDateSession{iSession};
        % dateSession = '20191113'; %setDateSession{iSession};
        
        dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
        
        fprintf(1, ':: analCa_ttestBPM.m :: processing %s_%s (#%d/%d) ...\n', nameSubj, dateSession, iSession, length(setDateSession));
        
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
        
        %% t-test to compare visual responses and baseline
        load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), 'tS_session_stim')
%         tS_session_stim = tS_session_stim(validIndCell, :);
        
        a = 0.1; %0.05;
        clear hval pval tval sd
        for iCell = 1:size(tS_session_stim, 1)
            for iStim = 1:size(tS_session_stim, 2)
                [h, p, ci, stats] = ttest(tS_session_stim(iCell, iStim).matAvgAmp_norm, tS_session_stim(iCell, iStim).matAvgAmp_b_norm, 'alpha', a/size(tS_session_stim,2));
                hval(iStim, iCell) = h;
                pval(iStim, iCell) = p;
                tval(iStim, iCell) = stats.tstat;
                sd(iStim, iCell) = stats.sd;
            end
        end
        
        paramTTest.alpha = a;
        paramTTest.alpha_BonferroniCorrection = a/size(tS_session_stim,2);
        paramTTest.pairedT_visualResp = 'tS_session_stim.matAvgAmp_norm';
        paramTTest.pairedT_baseline = 'tS_session_stim.matAvgAmp_b_norm';
        
        resultsTTest(iSession).dateSession = dateSession;
        resultsTTest(iSession).h = hval(:,validIndCell);
        resultsTTest(iSession).pval = pval(:,validIndCell);
        resultsTTest(iSession).tval = tval(:,validIndCell);
        resultsTTest(iSession).sd = sd(:,validIndCell);
        resultsTTest(iSession).flagVisualNeuron = sum(hval(:,validIndCell))>0;
        resultsTTest(iSession).indVisualNeuron = find(sum(hval(:,validIndCell))>0);
        
        if flagPlot            
            
            indNeuron = {}; %cell([]);
            indNeuron{1} = find(sum(resultsTTest(iSession).h)>5);
            indNeuron{2} = find(sum(resultsTTest(iSession).h)<6 & sum(resultsTTest(iSession).h)>1);
            indNeuron{3} = find(sum(resultsTTest(iSession).h)<2 & sum(resultsTTest(iSession).h)>0);
            indNeuron{4} = find(sum(resultsTTest(iSession).h)<1);
            cMap = jet(length(indNeuron));
            
            figure(100);clf;
            imagesc(neuron.Cn); colormap(gray);
            axis off
            hold on
            for iGroup = 1:length(indNeuron)
                plot(center(indNeuron{iGroup},2), center(indNeuron{iGroup},1), 'o', 'MarkerEdgeColor', cMap(iGroup,:), 'MarkerFaceColor', cMap(iGroup,:));
                hold on;
            end
            title(sprintf('%s %s: Paired T test results alpha0p1: broad (b) - intermediate (c) - very selective (y), not visual (r)', nameSubj, dateSession))
            print(gcf, fullfile(dirFig, sprintf('%s_%s_BPM_resultsPairedT_alpha0p1', nameSubj, dateSession)), '-depsc');
        end
        
    end
    
    save(fullfile(dirSave{iSubj}, 'BPM_pairedTTestResults_alpha0p1.mat'), 'paramTTest', 'resultsTTest')
    
end

% find(sum(hval)<1)


% Columns 1 through 15
%
%     11    12    13    14    16    17    18    20    22    29    32    38    44    52    53
%
%   Columns 16 through 26
%
%     58    70    72    77    83    85    98    99   107   109   116

%   Columns 1 through 15
%
%      7    13    14    17    18    32    33    38    46    58    71    72    83    85    99
%
%   Columns 16 through 17
%
%    107   109
