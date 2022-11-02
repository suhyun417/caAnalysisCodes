% saveDFL_sortedTimeSeries.m
%
% 2020/02/03 SHP

clear all; close all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/NIFVAULT/projects/parksh';
    dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/nifvault/projects/parksh';
    dirProcdata = '/nifvault/procdata/parksh';
    dirRawdata = '/nifvault/rawdata/parksh';
end

flagSave = 1;
flagPlot = 0;
flagSavePPTX = 0;

%% directory
setSubj = {'Tabla', 1; 'Max', 3};

for iSubj = 1:length(setSubj)
    
    nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
    FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    for iSession = 1:length(setDateSession)
        
        dateSession = setDateSession{iSession};
        % dateSession = '20191223';
        
        dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
        
        fprintf(1, 'Session #%d/%d (%s_%s): DFL runs are being processed....\n', iSession, length(setDateSession), dateSession, nameSubj);
        
        %% Read source data
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        
        %% Load the cell time series
        load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
        
        %% Stimulus timing info for DFL runs
        stimTiming_DFL = struct([]);
        for iRun = 1:length(tSeries_DFL)
            
            load(fullfile(dirProcdata_session, sprintf('DFL_%d_tML.mat', iRun)), 't_adj', 'stim', 'analog');
            
            %     stimTiming_DFL(iRun).filename_org = DataFile;
            stimTiming_DFL(iRun).t_adj = t_adj;
            
            stimTiming_DFL(iRun).stim = stim;
            %         stimTiming_DFL(iRun).indValidTrial = find(stim.trialError>1);
            
            fs = 10; %eventually should be retrieved from xml file instead of hard-coding
            %         locTrialStart = floor(t_adj.trialStart./1000./(1/fs));
            locMovieOn = floor(t_adj.movieOnset./(1000/fs))+1;
            locMovieOff = floor(t_adj.sendTTL_end./(1000/fs))+1;
            locReward = floor(t_adj.reward./(1000/fs))+1;
            %         locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./1000./(1/fs));
            %     t_onset_adj = t_onset/1000-delay;
            %     locStimOn = floor(t_onset_adj./(1/fs));
            %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
            
            %         stimTiming_DFL(iRun).locCaFrame.indTrial = find(stim.trialError>1);
            %         stimTiming_DFL(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_DFL(iRun).locCaFrame.indTrial);
            stimTiming_DFL(iRun).locCaFrame.durMovie = 120*fs;
            stimTiming_DFL(iRun).locCaFrame.locMovieOn = locMovieOn;
            stimTiming_DFL(iRun).locCaFrame.locMovieOff = locMovieOff;
            stimTiming_DFL(iRun).locCaFrame.locReward = locReward;
            
            stimTiming_DFL(iRun).analog = analog;
            
        end
        
        
        %% DFL: sort the timeseries for each cell and each movie
        catStim = cat(1, stimTiming_DFL(:).stim);
        catNameMovie = {catStim.nameMovie}';
        setMovie = unique(catNameMovie);
        
        tS_session = struct([]);
        for iMovie = 1:length(setMovie)
            setIndRun = find(contains(catNameMovie, setMovie{iMovie})>0);
            
            for iRun = 1:length(setIndRun)
                
                idRun = setIndRun(iRun);
                
                tS_session(iMovie).matTS_C_raw(:, :, iRun) = tSeries_DFL(idRun).C_raw(:, stimTiming_DFL(idRun).locCaFrame.locMovieOn:stimTiming_DFL(idRun).locCaFrame.locMovieOn+stimTiming_DFL(iRun).locCaFrame.durMovie)';
                tS_session(iMovie).matTS_C(:, :, iRun) = tSeries_DFL(idRun).C(:, stimTiming_DFL(idRun).locCaFrame.locMovieOn:stimTiming_DFL(idRun).locCaFrame.locMovieOn+stimTiming_DFL(iRun).locCaFrame.durMovie)';

                tS_session(iMovie).matTS_C_raw_zscore(:, :, iRun) = zscore(tS_session(iMovie).matTS_C_raw(:, :, iRun));
                tS_session(iMovie).avgTS_C_raw = median(tS_session(iMovie).matTS_C_raw(:, :, iRun), 3);
                tS_session(iMovie).avgTS_C_raw_zscore = median(tS_session(iMovie).matTS_C_raw_zscore(:, :, iRun), 3);
                
            end
            
            tS_session(iMovie).idStim = setMovie{iMovie};
            tS_session(iMovie).indRun = setIndRun;
        end

        %% save data into a mat file
        if flagSave
            % save file name
            saveFileName = fullfile(dirProcdata_session, 'DFL_ts_tML.mat');
            save(saveFileName, 'stimTiming_DFL', 'tS_*')
        end   
        
%         if flagPlot % 
%             figCheck = figure;
%             set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 880 315])
%             for iCell = 1:size(tS_session(1).matTS_C_raw, 2) %length(setCell) %1:size(matTS_movie, 1)
%                 figure(figCheck);
%                 SP(1) = subplot(2, 1, 1);
%                 imagesc(zscore(squeeze(matTS_movie(iCell, :, :)))')
%                 colormap(hot)
%                 ylabel(SP(1), 'Viewing')
%                 set(SP(1), 'CLim', [0 8]);
%                 %     plot(squeeze(matTS_movie(iCell, :, :)))
%                 SP(2) = subplot(2, 1, 2);
%                 plot(zscore(squeeze(matTS_movie(iCell, :, :))))
%                 axis tight
%                 title(SP(1), sprintf('Cell #%d/%d: movie %s', iCell, size(matTS_movie,1), tSeries_DFL(1), setMovie{iMovie}));
%                 
%                 set(SP(:), 'XTickLabel', 20:20:120)
%                 xlabel(SP(2), 'Time (s)');
%                 ylabel(SP(2), 'Norm. resp. (std)')
%                 
%                 F(iCell) = getframe(gcf);
%                 drawnow
%             end
%         end
%             
%         % save figures
%         if flagSavePPTX
%             fname_pptx = sprintf('%s_%s', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
%             exportFigsToPPTX(fname_pptx);
%             
%             switch lower(nameSubj)
%                 case 'tabla'
%                     dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
%                 case 'max'
%                     dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
%             end
%             movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
%         end
        
%         fname_movie = sprintf('%s_%s_DFL_%dCells_tSeries_movie%s', dateSession, nameSubj, length(setCell), setMovie{iMovie});
%         writerObj = VideoWriter(fullfile(dirProcdata_session, fname_movie), 'MPEG-4');
%         writerObj.FrameRate = 1;
%         open(writerObj);
%         for i = 1:length(F)
%             frame = F(i);
%             writeVideo(writerObj, frame);
%         end
%         close(writerObj);
    end
end
        
%         indNeuron = 1:size(matTS_movie, 1); %neuron.orderROIs('snr');
%         setNeuron = indNeuron;
%         figCheck = figure;
%         for iCell = 1:length(setNeuron) %1:size(matTS_movie, 1)
%             idCell = setNeuron(iCell);
%             figure(figCheck);
%             SP(1) = subplot(2, 1, 1);
%             imagesc(zscore(squeeze(matTS_movie(idCell, :, :)))')
%             colormap(hot)
%             set(SP(1), 'CLim', [0 8]);
%             %     plot(squeeze(matTS_movie(iCell, :, :)))
%             SP(2) = subplot(2, 1, 2);
%             plot(zscore(squeeze(matTS_movie(idCell, :, :))))
%             axis tight
%             title(SP(1), sprintf('Cell #%d/%d (Cell ID: %d): movie %s', iCell, size(matTS_movie,1), tSeries_DFL(1).idNeuron(idCell), setMovie{iMovie}));
%             
%             input('')
%         end
        
%         Reference
%         tS_session = struct([]);
%         tS_session(1).idRunTrial = idRunTrial;
%         tS_session(1).idStim = cat(1, stimTiming_BPM.idStim);
%         % tS_session(1).nameCondition = cat(1, stimTiming_BPM.infoStim);
%         tS_session(1).tS_trial = cat(2, tS_run.tS_trial); % Cell by Trial
%         
%         % Sort for different stimulus type
%         [sortStim, indTrialStim] = sort(tS_session.idStim);
%         setStim = unique(sortStim);
%         
%         tS_session_stim = struct([]);
%         for iCell = 1:size(tS_session(1).tS_trial, 1)
%             for iStim = 1:length(setStim)
%                 curStim = setStim(iStim);
%                 curIndTrial = indTrialStim(sortStim==curStim);
%                 
%                 tS_session_stim(iCell, iStim).idStim = curStim;
%                 tS_session_stim(iCell, iStim).indTrial = curIndTrial;
%                 tS_session_stim(iCell, iStim).indTrial_org = tS_session.idRunTrial(curIndTrial, :);
%                 tS_session_stim(iCell, iStim).matTS = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS);
%                 tS_session_stim(iCell, iStim).matTS_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS_norm);
%                 tS_session_stim(iCell, iStim).matAmp = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp);
%                 tS_session_stim(iCell, iStim).matAmp_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp_norm);
%                 tS_session_stim(iCell, iStim).matAvgAmp = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
%                 tS_session_stim(iCell, iStim).matAvgAmp_norm = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
%                 tS_session_stim(iCell, iStim).avgAmp = median(tS_session_stim(iCell, iStim).matAvgAmp); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
%                 tS_session_stim(iCell, iStim).avgAmp_norm = median(tS_session_stim(iCell, iStim).matAvgAmp_norm); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
%                 
%                 clear curStim curIndTrial
%             end
%         end
        