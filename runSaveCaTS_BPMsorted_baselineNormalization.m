% runSaveCaTS_BPMsorted_baselineNormalization
%
% 2025/11/13 created by SHP
% For Imaging Neuroscience revision
%   - Stimulus response from BPM session needs to be recomputed using
%   baseline subtraction & division by the baseline STD
%   - so instead of using currently saved normalized (simple z-scoring of the entire session) 
%   amplitude and baseline, re-compute and added it to the BPM_ts_tML.mat
%   - due to the save permission problem with KAIST linux machine & NAS
%   directory, I saved the results locally then manually copied the new
%   BPM_ts_tML.mat to procdata directory, and renamed the previous
%   BPM_ts_tML.mat to BPM_ts_tML_prev.mat
%   - after this code, runSaveCaTS_acrossSessionFOV_BPMsorted.m was run to
%   merge the matrix across sessions for each cell

clear all;



%% Directory settings
directory = setDir_shp;
dirProjects = directory.dirProjects;
dirProcdata = directory.dirProcdata;
dirRawdata = directory.dirRawdata;
dirFig = directory.dirFig;


addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));

%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; %1; %2; %1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%%
for iSession = 1:nSession

    clear tS* stim*

    dateSession = setDateSession{iSession}; %'20191205';

    dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');

    dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
    
    % load the time series
    load([dirProcdata, sprintf('/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)]); %, 'tS_session_stim', 'stimTiming_BPM')

    %% load each run and save normalized (basedline subtracted amplitude)
    for iRun = 1:length(tS_run)
        tS_trial_b_norm = struct([]);
        for iCell = 1:size(tS_run(iRun).tS_trial, 1)
            clear tS_trial cat_resp_b
            tS_trial = tS_run(iRun).tS_trial;
            cat_resp_b = cat(1, tS_trial(iCell, :).matAmp_b);
            b_mean = mean(cat_resp_b);
            b_std = std(cat_resp_b);

            % figure(100)
            % plot(cat_resp_b, 'o'); drawnow;
            % title(sprintf('Run #%d, Cell #%02d: mean %2.02f, std %2.02f ', iRun, iCell, b_mean, b_std))
            % input('')

            for iTrial = 1:size(tS_trial, 2)
                clear mat*
                % normalize with baseline (basedline mean subtraction & divided by baseline std)
                matTS_org = tS_trial(iCell, iTrial).matTS;
                matTS_b_norm = (matTS_org - b_mean)./b_std;
                matAmp_b_norm = matTS_b_norm(tS_trial(iCell, iTrial).locRespAmp);

                tS_trial_b_norm(iCell, iTrial).matTS_b_norm = matTS_b_norm;
                tS_trial_b_norm(iCell, iTrial).b_mean = b_mean;
                tS_trial_b_norm(iCell, iTrial).b_std = b_std;
                tS_trial_b_norm(iCell, iTrial).locRespAmp = tS_trial(iCell, iTrial).locRespAmp; %stimulus response window
                tS_trial_b_norm(iCell, iTrial).resp_b_norm = matAmp_b_norm; %stimulus response
                tS_trial_b_norm(iCell, iTrial).resp_b_norm_avg = mean(matAmp_b_norm); %stimulus response

                %     figure(150); cla;
                %     plot(matTS_org); hold on
                %     plot(matTS_b_norm, 'r-'); % after normalization
                %     title(sprintf('Run #%d, Cell #%02d: mean %2.02f, std %2.02f, stim %d', iRun, iCell, b_mean, b_std, tS_trial(iCell, iTrial).idStim))
                %     input('')

            end
        end
        tS_run(iRun).tS_trial_b_norm = tS_trial_b_norm;
    end
    clear tS_trial*

    % merge across runs
    % merge info for session
    % % record the Run identity for later
    % idRunTrial = [];
    % for iRun = 1:length(stimTiming_BPM)
    %     idRunTrial = cat(1, idRunTrial, cat(2, ones(size(stimTiming_BPM(iRun).indValidTrial)).*iRun, [1:length(stimTiming_BPM(iRun).indValidTrial)]', ...
    %         stimTiming_BPM(iRun).indValidTrial, stimTiming_BPM(iRun).indValidTrial_orgBhv));
    % end
    % tS_session = struct([]);
    % tS_session(1).idRunTrial = idRunTrial;
    % tS_session(1).idRunTrial_description = 'Run number, Trial index in this file within each Run, Trial index in tML.mat file, Trial index in original bhv file';
    % tS_session(1).idStim = cat(1, stimTiming_BPM.idStim);
    % tS_session(1).nameCondition = cat(1, stimTiming_BPM.infoStim);
    tS_session(1).tS_trial_b_norm = cat(2, tS_run.tS_trial_b_norm); % Cell by Trial

    % Sort for different stimulus type
    [sortStim, indTrialStim] = sort(tS_session.idStim);
    setStim = unique(sortStim);

    for iCell = 1:size(tS_session(1).tS_trial_b_norm, 1)
        for iStim = 1:length(setStim)
            curStim = setStim(iStim);
            curIndTrial = indTrialStim(sortStim==curStim);

            %         tS_session_stim(iCell, iStim).idStim = curStim;
            %         tS_session_stim(iCell, iStim).indTrial = curIndTrial;
            %         tS_session_stim(iCell, iStim).indTrial_org = tS_session.idRunTrial(curIndTrial, :);

            % baseline-normalized values
            tS_session_stim(iCell, iStim).baselineNorm_param = '20251118_baseline from the entire run is used to normalize and termed "baselineNorm"';
            tS_session_stim(iCell, iStim).matTS_baselineNorm = cat(2, tS_session.tS_trial_b_norm(iCell, curIndTrial).matTS_b_norm); % entire trial TS
            tS_session_stim(iCell, iStim).matResp_baselineNorm = cat(2, tS_session.tS_trial_b_norm(iCell, curIndTrial).resp_b_norm); % stimulus evoked responses
            tS_session_stim(iCell, iStim).matRespAvg_baselineNorm = cat(1, tS_session.tS_trial_b_norm(iCell, curIndTrial).resp_b_norm_avg); % response amplitudes

            %         tS_session_stim(iCell, iStim).matTS = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS);
            %         %                 tS_session_stim(iCell, iStim).matTS_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matTS_norm);
            %         % amplitude
            %         tS_session_stim(iCell, iStim).matAmp = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp);
            %         %                 tS_session_stim(iCell, iStim).matAmp_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp_norm);
            %         tS_session_stim(iCell, iStim).matAvgAmp = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
            %         %                 tS_session_stim(iCell, iStim).matAvgAmp_norm = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
            %         tS_session_stim(iCell, iStim).avgAmp = median(tS_session_stim(iCell, iStim).matAvgAmp); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
            %         %                 tS_session_stim(iCell, iStim).avgAmp_norm = median(tS_session_stim(iCell, iStim).matAvgAmp_norm); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);
            %         % baseline
            %         tS_session_stim(iCell, iStim).matAmp_b = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp_b);
            %         %                 tS_session_stim(iCell, iStim).matAmp_b_norm = cat(2, tS_session.tS_trial(iCell, curIndTrial).matAmp_b_norm);
            %         tS_session_stim(iCell, iStim).matAvgAmp_b = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_b);
            %         tS_session_stim(iCell, iStim).matAvgAmp_b_norm = cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_b_norm);
            %         tS_session_stim(iCell, iStim).avgAmp_b = median(tS_session_stim(iCell, iStim).matAvgAmp_b); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp);
            %         %                 tS_session_stim(iCell, iStim).avgAmp_b_norm = median(tS_session_stim(iCell, iStim).matAvgAmp_b_norm); %cat(1, tS_session.tS_trial(iCell, curIndTrial).avgAmp_norm);

            clear curStim curIndTrial
        end
    end

    % save file name
    saveFileName = sprintf('/home/parks23/Research/0Marmoset/Ca/tempData/%s_%s_BPM_ts_tML.mat', nameSubj, dateSession);
%     saveFileName = fullfile(dirProcdata_session, 'BPM_ts_tML.mat');
    save(saveFileName, 'stimTiming_BPM', 'tS*')

end

