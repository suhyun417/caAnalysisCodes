% analCa_BPMDFL_allRunsSourceExtract.m
%
% 2020/01/28 SHP

clear all; close all;

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

%% directory
setNameSubj = {'Tabla', 'Max'};

for iSubj = 1:length(setNameSubj)
    
    nameSubj = setNameSubj{iSubj}; %'Tabla'; %'Max';
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    for iSession = 1:length(setDateSession)
        
        close all;
        
        dateSession = setDateSession{iSession};
        % dateSession = '20191223';
        
        
        dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
        
        
        
        %% Read source data
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        
        %% Load the cell time series
        load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
        load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
        
        
        %% Stimulus timing info for BPM runs
        stimTiming_BPM = struct([]);
        for iRun = 1:length(tSeries_BPM)
            
            load(fullfile(dirProcdata_session, sprintf('BPM_%d_tML.mat', iRun)), 't_adj', 'infoTrial', 'analog')
            
            validT = find(t_adj.trialStart>0);
            
            %     stimTiming_BPM(iRun).t_adj = t_adj;
            stimTiming_BPM(iRun).condMat = infoTrial.condMat(validT, :);
            stimTiming_BPM(iRun).idStim = infoTrial.idStim(validT);
            stimTiming_BPM(iRun).setIdStim = infoTrial.uniqueSetIdStim;
            stimTiming_BPM(iRun).infoStim = infoTrial.infoStim';
            
            fs = 10; %eventually should be retrieved from xml file instead of hard-coding
            locFrameTrialStart = floor(t_adj.trialStart(validT)./(1000/fs))+1; % add one because this is frame index and the first frame is "1"
            locFrameStimOn = floor(t_adj.stimOnset(validT)./(1000/fs))+1;
            locFrameBlankOn_afterStim = floor(t_adj.blankOnset_afterStim(validT)./(1000/fs))+1;
            
            stimTiming_BPM(iRun).indValidTrial = validT;
            stimTiming_BPM(iRun).locCaFrame.trialStart = locFrameTrialStart;
            stimTiming_BPM(iRun).locCaFrame.stimOn = locFrameStimOn;
            stimTiming_BPM(iRun).locCaFrame.blankOn_afterStim = locFrameBlankOn_afterStim;
            
            stimTiming_BPM(iRun).analog = analog;
            
        end
        
        %% BPM: sort the timeseries for each cell and each condition
        tS_run = struct([]);
        for iRun = 1:size(tSeries_BPM, 2)
            
            matTS = tSeries_BPM(iRun).C_raw';
            matTS_norm = zscore(matTS);
            
            tS_trial = struct([]);
            for iCell = 1:size(matTS, 2)
                for iTrial = 1:length(stimTiming_BPM(iRun).indValidTrial)
                    
                    tS_trial(iCell, iTrial).matTS = matTS(stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)-10:stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+35, iCell); %
                    tS_trial(iCell, iTrial).matTS_norm = matTS_norm(stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)-10:stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+35, iCell);
                    
                    win_ms = [500 1500]; % time window of amplitude calculation (in ms: zero is stim Onset)
                    fs = 10;
                    win_frame = win_ms./(1000/fs);
                    
                    tS_trial(iCell, iTrial).idStim = stimTiming_BPM(iRun).idStim(iTrial);
                    tS_trial(iCell, iTrial).amp_win_ms = win_ms;
                    tS_trial(iCell, iTrial).amp_win_frame = win_frame;
                    
                    tS_trial(iCell, iTrial).matAmp = matTS(stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+win_frame(1): stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+win_frame(2), iCell);
                    tS_trial(iCell, iTrial).avgAmp = median(tS_trial(iCell, iTrial).matAmp);
                    tS_trial(iCell, iTrial).matAmp_norm = matTS_norm(stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+win_frame(1): stimTiming_BPM(iRun).locCaFrame.stimOn(iTrial)+win_frame(2), iCell);
                    tS_trial(iCell, iTrial).avgAmp_norm = median(tS_trial(iCell, iTrial).matAmp_norm);
                    
                end
            end
            
            tS_run(iRun).tS_trial = tS_trial;
        end
        
        % merge across runs
        % merge info for session
        % record the Run identity for later
        idRunTrial = [];
        for iRun = 1:length(stimTiming_BPM)
            idRunTrial = cat(1, idRunTrial, cat(2, ones(size(stimTiming_BPM(iRun).indValidTrial)).*iRun, stimTiming_BPM(iRun).indValidTrial));
        end
        tS_session = struct([]);
        tS_session(1).idRunTrial = idRunTrial;
        tS_session(1).idStim = cat(1, stimTiming_BPM.idStim);
        % tS_session(1).nameCondition = cat(1, stimTiming_BPM.infoStim);
        tS_session(1).tS_trial = cat(2, tS_run.tS_trial); % Cell by Trial
        
        % Sort for different stimulus type
        [sortStim, indTrialStim] = sort(tS_session.idStim);
        setStim = unique(sortStim);
        
        tS_session_stim = struct([]);
        for iCell = 1:size(tS_session(1).tS_trial, 1)
            for iStim = 1:length(setStim)
                curStim = setStim(iStim);
                curIndTrial = indTrialStim(sortStim==curStim);
                
                tS_session_stim(iCell, iStim).idStim = curStim;
                tS_session_stim(iCell, iStim).indTrial = curIndTrial;
                tS_session_stim(iCell, iStim).indTrial_org = tS_session.idRunTrial(curIndTrial, :);
                tS_session_stim(iCell, iStim).matTS = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS);
                tS_session_stim(iCell, iStim).matTS_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS_norm);
                tS_session_stim(iCell, iStim).matAmp = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp);
                tS_session_stim(iCell, iStim).matAmp_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp_norm);
                tS_session_stim(iCell, iStim).matAvgAmp = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
                tS_session_stim(iCell, iStim).matAvgAmp_norm = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
                tS_session_stim(iCell, iStim).avgAmp = median(tS_session_stim(iCell, iStim).matAvgAmp); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
                tS_session_stim(iCell, iStim).avgAmp_norm = median(tS_session_stim(iCell, iStim).matAvgAmp_norm); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
                
                clear curStim curIndTrial
            end
        end
        
        % save file name
        saveFileName = fullfile(dirProcdata_session, 'BPM_ts_tML.mat');
        save(saveFileName, 'stimTiming_BPM', 'tS*')
        
        % % across run merge
        % tempCat = cat(3, tS_run.tS_stim); % Cell by Stimulus by Run
        % tS_session = struct([]);
        % for iCell = 1:size(tempCat,1)
        %     for iStim = 1:size(tempCat, 2)
        %
        %         tS_session(iCell, iStim).idStim = tempCat(iCell, iStim, 1).idStim;
        %         tS_session(iCell, iStim).matTS = cat(2, tempCat(iCell, iStim, :).matTS);
        %         tS_session(iCell, iStim).matTS_norm = cat(2, tempCat(iCell, iStim, :).matTS_norm);
        %         tS_session(iCell, iStim).matAmp = cat(2, tempCat(iCell, iStim, :).matAmp);
        %         tS_session(iCell, iStim).matAmp_norm = cat(2, tempCat(iCell, iStim, :).matTS_norm);
        %         tS_session(iCell, iStim).matAvgAmp = cat(1, tempCat(iCell, iStim, :).avgAmp);
        %         tS_session(iCell, iStim).matAvgAmp_norm = cat(1, tempCat(iCell, iStim, :).avgAmp_norm);
        %         tS_session(iCell, iStim).avgAmp = median(tS_session(iCell, iStim).matAvgAmp);
        %         tS_session(iCell, iStim).avgAmp_norm = median(tS_session(iCell, iStim).matAvgAmp_norm);
        %
        %         % record the Run identity for later
        %         idRunTrial = [];
        %         for iRun = 1:size(tempCat, 3)
        %             idRunTrial = cat(1, idRunTrial, cat(2, ones(size(tempCat(iCell, iStim, iRun).indTrial)).*iRun, tempCat(iCell, iStim, iRun).indTrial));
        %         end
        %         tS_session(iCell, iStim).idRunTrial = idRunTrial;
        %
        %     end
        % end
        
        %% FUN TIME
        %
        
        condName_BPM = {'human face', 'marmoset face', 	'marmoset body', 'scene', 'non familiar object', 'hands and catcher',...
            'phase scrambled', 'space scrambled', 'grating', 'random dot motion'};
        
        setCond = [1 2 5 6 10]; %cat(2, 11:15, 21:25, 51
        nImage = 5;
        % temp = repmat([1:nImage], length(setCond), 1);
        % tempCondMat = cat(2, repmat(setCond', nImage, 1), temp(:));
        % setStimID = sort(tempCondMat(:,1)*10 + tempCondMat(:,2));
        setCondName = condName_BPM(setCond); %{infoTrial.infoStim([1:6:25]).nameCondition};
        
        neuron_b = neuron.batches{1}.neuron;
        thr = 0.3; % the lower the smaller (more centralized) the contour
        Coor = neuron_b.get_contours(thr);
        imgFOV = neuron_b.Cn.*neuron_b.PNR;
        
%         figure;
        neuron_b.show_contours([], [], imgFOV, 'true');
        
        
        % summary responses
        matAmpCellStim = reshape(cat(1, tS_session_stim.avgAmp), size(tS_session_stim));
        figure;
        set(gcf, 'Color', 'w', 'Position', [680         602        1080         376])
        imagesc(matAmpCellStim')
        colormap(jet)
        set(gca, 'CLim', [-1 1].*2)
        xlabel('Cells')
        ylabel('Stimulus')
%         condName = {infoTrial.infoStim([1:6:25]).nameCondition};
        set(gca, 'YTickLabel', setCondName)
        
        % fig_BPM = figure;
        % set(fig_BPM, 'Color', 'w', 'Position', [67          58        1144         917]);
        cMap_run = bone(length(tSeries_BPM)+1); % for each run
        cMap_cond = cool(length(setCond));
        catStimIDSession = cat(1, tS_session_stim(1,:).idStim);
        for iCell = 1:size(tS_session_stim,1)
            %     iCell = orderCell(iC);
            %     figure(fig_BPM); clf;
            fig_BPM = figure;
            set(fig_BPM, 'Color', 'w', 'Position', [67          58        1144         917]);
            %     ylim = [];
            clear sp_stim
            for iCond = 1:length(setCond)
                idCond = setCond(iCond);
                for iM = 1:nImage
                    indStim = find(catStimIDSession == 10*idCond + iM);
                    
                    figure(fig_BPM);
                    sp_stim(iCond, iM) = subplot(length(setCond)+1, nImage+1, (iCond-1)*(nImage+1)+iM);
                    set(gca, 'ColorOrder', cMap_run(tS_session_stim(iCell, indStim).indTrial_org(:,1), :));
                    hold on;
                    plot([-1:0.1:3.5], tS_session_stim(iCell, indStim).matTS_norm, 'LineWidth', 1);
                    set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
                    
                    avgTS(:,iM) = mean(tS_session_stim(iCell, indStim).matTS_norm, 2);
                    
                    if iM==1
                        ylabel(sp_stim(iCond, iM), setCondName{iCond}, 'Color', cMap_cond(iCond, :));
                    end
                end
                
                figure(fig_BPM);
                sp_stim(iCond, nImage+1) = subplot(length(setCond)+1, nImage+1, (iCond-1)*(nImage+1)+nImage+1);
                cla;
                set(gca, 'ColorOrder', jet(nImage).*0.8);
                hold on;
                plot([-1:0.1:3.5], avgTS)
                set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
                
                avgTS_cond(:, iCond) = mean(avgTS, 2);
            end
            
            sp_stim(length(setCond)+1, 1) = subplot(length(setCond)+1, nImage+1, length(setCond)*(nImage+1)+1);
            imagesc(imgFOV)
            colormap(sp_stim(length(setCond)+1, 1), gray)
            hold on
            plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'c.')
            axis off
            
            sp_stim(length(setCond)+1, nImage+1) = subplot(length(setCond)+1, nImage+1, (length(setCond)+1)*(nImage+1));
            set(gca, 'ColorOrder', cMap_cond);
            hold on;
            plot([-1:0.1:3.5], avgTS_cond)
            set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
            ylabel('norm. resp (z)')
            xlabel('Time from stim on (s)')
            
            title(sp_stim(1, 3), sprintf('Cell #%d/%d', iCell, size(tS_session_stim,1)))
            
            xlabel(sp_stim(5, 3), 'Time from stim on (s)')
            
            %     input('')
        end
        
        % save figures 
        fname_pptx = sprintf('%s_%s', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
        exportFigsToPPTX(fname_pptx);
        
        switch lower(nameSubj)
            case 'tabla'
                dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
            case 'max'
                dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
        end
        movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
        
    end
end

%
%
%     for iStim = 1:length(setStimID)
%
%         indStim = find(catStimIDSession == setStimID(iStim));
%
%         figure(fig_BPM);
%         SP(iStim) = subplot(length(setCond)+1, nImage, iStim);
%
%         % each trial in heatmap
%         imagesc(tS_session_stim(iCell, indStim).matTS_norm')
%         colormap(jet)
%         set(gca, 'CLim', [-1 1].*3)
%
%         avgTS(:,iStim) = median(tS_session_stim(iCell, iStim).matTS_norm, 2);
%
% %         % mean across trial
% %         hold on;
% %         plot([-1:0.1:3], median(tS_session_stim(iCell, iStim).matTS_norm, 2), 'k-', 'LineWidth', 2);
% %         set(gca, 'XLim', [-1 4], 'XTick', -1:1:4)
%
% %         % each trial
% %         set(gca, 'ColorOrder', cMap(tS_session_stim(iCell, indStim).indTrial_org(:,1), :));
% %         hold on;
% %         plot([-1:0.1:3], tS_session_stim(iCell, indStim).matTS_norm, 'LineWidth', 1);
% %         set(gca, 'XLim', [-1 3], 'XTick', -1:1:3)
%
% %         ylim = [-2 4]; %cat(1, ylim, get(SP(iStim), 'YLim'));
%     end
%
% %     set(SP, 'YLim', [min(ylim(:,1)) max(ylim(:,2))])
%     title(SP(3), sprintf('Cell #%d/%d', iCell, size(tS_session_stim,1)))
%     ylabel(SP(1), setCondName{1}); ylabel(SP(6), setCondName{2}); ylabel(SP(11), setCondName{3});
%     ylabel(SP(16), setCondName{4}); ylabel(SP(21), setCondName{5});
%     xlabel(SP(23), 'Time from stim on (s)')
%
%     SP((length(setCond)*nImage)+1) = subplot(length(setCond)+1, nImage, (length(setCond)*nImage)+1);
%     imagesc(imgFOV)
%     colormap(SP((length(setCond)*nImage)+1), gray)
%     hold on
%     plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'c.')
%     axis off
%
%     input('')
% end


% figure
% imagesc(neuron_b.Cn.*neuron_b.PNR)
% colormap(gray)
% hold on
% iCell = 1;
% plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'co')
%
% matAmpCellStim = reshape(cat(1, tS_session_stim.avgAmp), size(tS_session_stim));
% figure;
% set(gca, 'Color', 'w', 'Position', [680         602        1080         376])
% imagesc(matAmpCellStim')
% colormap(jet)
% set(gca, 'CLim', [-1 1].*2)
% xlabel('Cells')
% ylabel('Stimulus')
% condName = {infoTrial.infoStim([1:6:25]).nameCondition};
% set(gca, 'YTickLabel', condName)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bis hier
% %%
% fig_selectivity = figure;
%
% for iCell = 1:size(tSeries_BPM_Cond, 1)
%     fprintf(1, 'Cell #%d: ', iCell)
%
%     %     iRun = 1;
%     for iRun = 1:size(tSeries_BPM_Cond,2)
%         %         curTS_catCond = [];
%         %         curTS_catCond = cat(2, tSeries_BPM_Cond(iCell, iRun, :).matTS);
%
%         for iCond = 1:5
%             figure(fig_selectivity)
%             SP(5*(iRun-1)+iCond) = subplot(size(tSeries_BPM_Cond,2), 5, 5*(iRun-1)+iCond);
%             plot(tSeries_BPM_Cond(iCell, iRun, iCond).matTS_norm)
%             %             plot(curTS_catCond(:,:,iCond))
%             ylim([-3 10])
%         end
%
%         title(SP(3), sprintf('Cell #%d', iCell))
%
%     end
%
%     input('')
% end
%
% for iRun = 1:size(tSeries_BPM, 2)
%
%     %     indTrial = stimTiming_BPM(iRun).locCaFrame.indTrial;
%     %     indValidTrial = find(stimTiming_BPM(iRun).locCaFrame.locTrialStart>0);
%     %     condMat = stimTiming_BPM(iRun).stim.condMat(stimTiming_BPM(iRun).locCaFrame.indTrial(indValidTrial),:);
%
%     [sortedCond, indTrialCond] = sortrows(stimTiming_BPM(iRun).condMat);
%     catMatAmp = cat(2, resultsTrial(indTrialCond).matAmp);
%     catMnAmp = cat(2, resultsTrial(indTrialCond).mnAmp);
%     % figure;
%     % plot(catMatAmp(5,:), 'o');
%     % figure;
%     % plot(catMatAmp(14,:), 'o');  hold on
%     figure(100);
%     for iC = 1:size(catMatAmp, 1)
%         figure(100);cla;
%         plot(catMatAmp(iC,:), 'o');  hold on
%         line([110 209 319 429; 110 209 319 429], [-3 -3 -3 -3; 7 7 7 7], 'Color', 'm')
%         title(sprintf('Cell #%d/%d:', iC, size(catMatAmp, 1)))
%         input('')
%     end
% end
%
% clear tempCat
% setCell = 1:size(tSeries_BPM_Cond, 1); %[5 14 37 38 46 47 48 54 76 77]; %1:size(tSeries_BPM_Cond, 1); %[5 14 37 38 46 47 48]; % [1 4 6];
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell); %4; %1;
%     for iCond = 1:5
%         tempCat{iCond} = cat(2, tSeries_BPM_Cond(idCell, :, iCond).matTS);
%         tempCat_norm{iCond} = cat(2, tSeries_BPM_Cond(idCell, :, iCond).matTS_norm);
%     end
%     for iCond = 1:5
%         mnTS(:,iCond) = mean(tempCat{iCond}, 2);
%         steTS(:,iCond) = std(tempCat{iCond}, [], 2)./sqrt(size(tempCat{iCond}, 2));
%
%         mnTS_norm(:,iCond) = mean(tempCat_norm{iCond}, 2);
%     end
%
%     figure(100);
%     %     set(gcf, 'Position', [100 100 800 500])
%     %     SP(1) = subplot(1,2,1);
%     %     plot([1:51].*0.1, mnTS, 'LineWidth', 2)
%     %     legend('HF', 'MF', 'NO', 'FO', 'RD')
%     %     title(sprintf('Cell #%d', iCell))
%     %     axis tight
%
%     %     SP(2) = subplot(1,2,2);
%     plot([-1:0.1:4], mnTS_norm, 'LineWidth', 2)
%     legend('HF', 'MF', 'NO', 'FO', 'RD')
%     title(sprintf('Cell #%d', iCell))
%     ylabel('zscore')
%     set(gca, 'YLim', [-1 4], 'XLim', [-1 4], 'XTick', -1:1:4)
%     %     axis tight
%
%     set(gca, 'Box', 'off', 'TickDIr', 'out', 'LineWidth', 2)
%
%     input('')
%     %     print(gcf, fullfile(dirFig, sprintf('%s_%s_ConditionMerge_Cell%d', dateSession, nameSubj, iCell)), '-depsc')
%
% end
%
% % Draw a contour (this is more manual way than neuron.show_contours();)
% thr = 0.12; % the lower the smaller (more centralized) the contour
% Coor = neuron_b.get_contours(thr);
% figure
% imagesc(neuron_b.Cn.*neuron_b.PNR)
% colormap(gray)
% hold on
% iCell = 1;
% plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'co')
%
%
% %% Stimulus timing info for DFL runs
% for iRun = 1:length(tSeries_DFL)
%
%     load(fullfile(dirProcdata_session, sprintf('DFL_%d_tML.mat', iRun)), 't_adj', 'stim');
%
%     %     stimTiming_DFL(iRun).filename_org = DataFile;
%     stimTiming_DFL(iRun).t_adj = t_adj;
%
%     stimTiming_DFL(iRun).stim = stim;
%     %         stimTiming_DFL(iRun).indValidTrial = find(stim.trialError>1);
%
%     fs = 10; %eventually should be retrieved from xml file instead of hard-coding
%     %         locTrialStart = floor(t_adj.trialStart./1000./(1/fs));
%     locMovieOn = floor(t_adj.movieOnset./(1000/fs))+1;
%     locMovieOff = floor(t_adj.sendTTL_end./(1000/fs))+1;
%     locReward = floor(t_adj.reward./(1000/fs))+1;
%     %         locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./1000./(1/fs));
%     %     t_onset_adj = t_onset/1000-delay;
%     %     locStimOn = floor(t_onset_adj./(1/fs));
%     %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
%
%     %         stimTiming_DFL(iRun).locCaFrame.indTrial = find(stim.trialError>1);
%     %         stimTiming_DFL(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_DFL(iRun).locCaFrame.indTrial);
%     stimTiming_DFL(iRun).locCaFrame.locMovieOn = locMovieOn;
%     stimTiming_DFL(iRun).locCaFrame.locMovieOff = locMovieOff;
%     stimTiming_DFL(iRun).locCaFrame.locReward = locReward;
%
% end
%
%
% %% DFL: sort the timeseries for each cell and each movie
% catStim = cat(1, stimTiming_DFL(:).stim);
% catNameMovie = {catStim.nameMovie}';
% setMovie = unique(catNameMovie);
%
% iMovie = 1;
% setIndRun = find(contains(catNameMovie, setMovie{iMovie})>0);
%
% clear matTS_movie
% for iRun = 1:length(setIndRun)
%
%     idRun = setIndRun(iRun);
%     matTS_movie(:, :, iRun) = tSeries_DFL(idRun).C_raw(:, stimTiming_DFL(idRun).locCaFrame.locMovieOn:stimTiming_DFL(idRun).locCaFrame.locMovieOff);
% end
%
% indNeuron = 1:size(matTS_movie, 1); %neuron.orderROIs('snr');
% setNeuron = indNeuron;
% figCheck = figure;
% for iCell = 1:length(setNeuron) %1:size(matTS_movie, 1)
%     idCell = setNeuron(iCell);
%     figure(figCheck);
%     SP(1) = subplot(2, 1, 1);
%     imagesc(zscore(squeeze(matTS_movie(idCell, :, :)))')
%     colormap(hot)
%     set(SP(1), 'CLim', [0 8]);
%     %     plot(squeeze(matTS_movie(iCell, :, :)))
%     SP(2) = subplot(2, 1, 2);
%     plot(zscore(squeeze(matTS_movie(idCell, :, :))))
%     axis tight
%     title(SP(1), sprintf('Cell #%d/%d (Cell ID: %d): movie %s', iCell, size(matTS_movie,1), tSeries_DFL(1).idNeuron(idCell), setMovie{iMovie}));
%
%     input('')
% end