% temp_analyzeDFL.m
%
% quick & dirty analysis of DFL sessions
% just checking the data and try out bunch of stuff, before making the code
% neat and tidy

clear all; close all;

%% directory
nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);
dirPreproc = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession, '_preproc');

dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %

% list of DFL imaging files
listRun_DFL = {'131229', '131728', '132314', '132823', '133237', '133624'}; %20191113_Tabla

%% Copy DFF data from /_preproc  with renaming of "DFL_1" 
d_file_imaging = dir(fullfile(dirProcdata_session, 'DFL*_dFF.mat'));
if isempty(d_file_imaging)
    d_dFF = dir(fullfile(dirPreproc, '*_dFF.mat'));
    listDFF = {d_dFF.name}';
    for iRun = 1:length(listRun_DFL)
        indFile = contains(listDFF, listRun_DFL{iRun});
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(dirPreproc, listDFF{indFile}), fullfile(dirProcdata_session, sprintf('DFL_%d_%s_dFF.mat', iRun, listRun_DFL{iRun})))
    end
else
    fprintf(' *_dFF.mat files are already present in %s\n', dirProcdata_session)
end

%% Read source data
addpath('/projects/parksh/_toolbox/CNMF_E/'); 
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_DFL*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); % manually copied and renamed the file this time
% Deleted 5 neurons:  2, 48, 74, 75, 78, from Tabla's Sources2D file
% (Sources2D_02-Dec_10_50_16.mat)

% indOrder_neuron = neuron.orderROIs('snr');
% matSource_org = neuron.reshape(neuron.A, 2); % 3D structure of spatial information of neurons
% nNeuron_org = size(matSource_org, 3);


%% Read each run's DFF and timing files
if ~exist(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'), 'file')
    
    d_file_imaging = dir(fullfile(dirProcdata_session, 'DFL*_dFF.mat'));
%     d_file_timing = dir(fullfile(dirProcdata_session, 'DFL*_tML.mat'));
    
% % Getting tSeries from dFF files computed outside of CNMFe
%     for iRun = 1:length(listRun_DFL)
%         
%         load(fullfile(dirProcdata_session, d_file_imaging(iRun).name))
%         % matTS = NaN(size(Y_dFF, 3), size(matSource_org,3)); % time x cell
%         matDFF = reshape(Y_dFF.*100, size(Y_dFF,1)*size(Y_dFF,2), size(Y_dFF,3)); % change it to percent and convert 3D to 2D        
% 
%         for iCell = 1:size(neuron.A, 2)
%             
%             clear tempTS
%             tempTS = matDFF(neuron.A(:,iCell)>0, :);
%             
%             tSeries_DFL(iCell, iRun).idNeuron_org =neuron.ids(iCell);
%             tSeries_DFL(iCell, iRun).matTS = tempTS;
%             tSeries_DFL(iCell, iRun).nPixel = size(tempTS, 1);
%             tSeries_DFL(iCell, iRun).mnTS = mean(tempTS);
%             tSeries_DFL(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
%             
%         end        
%         clear Y_dFF matDFF
%         
%     end
    
    % Using adjusted (ring background subtracted) time series from CNMFe
    count = 0;
    for iRun = 1:length(listRun_DFL)
        
        % get the # of frames from each session without loading the matrix
        matObj = matfile(fullfile(dirProcdata_session, d_file_imaging(iRun).name));
        [d1 d2 nFrame] = size(matObj, 'Y');
        
        tSeries_DFL(iRun).idNeuron_org = neuron.ids;
        tSeries_DFL(iRun).mnTS_C_raw = neuron.C_raw(:, count+1:count+nFrame);  %mean(tempTS);
        tSeries_DFL(iRun).mnTS_C = neuron.C(:, count+1:count+nFrame);
        
        count = count + nFrame;


%         for iCell = 1:size(neuron.A, 2)
% %             clear tempTS
% %             tempTS = matDFF(neuron.A(:,iCell)>0, :);
%             
%             tSeries_DFL(iCell, iRun).idNeuron_org =neuron.ids(iCell);
% %             tSeries_DFL(iCell, iRun).matTS = tempTS;
% %             tSeries_DFL(iCell, iRun).nPixel = size(tempTS, 1);
%             tSeries_DFL(iCell, iRun).mnTS_C_raw = neuron.C_raw(iCell, count+1:count+nFrame);  %mean(tempTS);
%             tSeries_DFL(iCell, iRun).mnTS_C = neuron.C(iCell, count+1:count+nFrame);
% %             tSeries_DFL(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
%             
%         end   
%         count = count+nFrame;
    end
        
   resultsCNMFe.A = neuron.A;
   resultsCNMFe.options = neuron.options;
   resultsCNMFe.Cn = neuron.Cn;
   resultsCNMFe.PNR = neuron.PNR;
   resultsCNMFe.ids = neuron.ids;
%    resultsCNMFe.indOrderNeuron = indOrder_neuron;
    
    d_file_timing = dir(fullfile(dirProcdata_session, 'DFL*_tML.mat'));
    
    for iRun = 1:length(listRun_DFL)
        
        load(fullfile(dirProcdata_session, d_file_timing(iRun).name), 'DataFile', 't_adj', 'stim')
        
        stimTiming_DFL(iRun).filename_org = DataFile;
        stimTiming_DFL(iRun).t_adj = t_adj;
        
        stimTiming_DFL(iRun).stim = stim;
%         stimTiming_DFL(iRun).indValidTrial = find(stim.trialError>1);
        
        fs = 10; %eventually should be retrieved from xml file instead of hard-coding
%         locTrialStart = floor(t_adj.trialStart./1000./(1/fs));
        locMovieOn = floor(t_adj.movieOnset./1000./(1/fs));
        locMovieOff = floor(t_adj.sendTTL_end./1000./(1/fs));
        locReward = floor(t_adj.reward./1000./(1/fs));
%         locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./1000./(1/fs));
        %     t_onset_adj = t_onset/1000-delay;
        %     locStimOn = floor(t_onset_adj./(1/fs));
        %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
        
%         stimTiming_DFL(iRun).locCaFrame.indTrial = find(stim.trialError>1);
%         stimTiming_DFL(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_DFL(iRun).locCaFrame.indTrial);
        stimTiming_DFL(iRun).locCaFrame.locMovieOn = locMovieOn;
        stimTiming_DFL(iRun).locCaFrame.locMovieOff = locMovieOff;
        stimTiming_DFL(iRun).locCaFrame.locReward = locReward;
        
        
    end
    
    % interim save
    save(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'), 'tSeries_DFL', 'stimTiming_DFL', 'resultsCNMFe')
    
else % if the file already exists
    load(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'))
    fprintf(1, 'Time series data in %s is loaded to the workspace. \n', fullfile(dirProcdata_session, 'DFL_ts_tML.mat'))
end


%% sort the timeseries for each cell and each movie
catStim = cat(1, stimTiming_DFL(:).stim);
catNameMovie = {catStim.nameMovie}';
setMovie = unique(catNameMovie);

iMovie = 1;
setIndRun = find(contains(catNameMovie, setMovie{iMovie})>0);

clear matTS_movie
for iRun = 1:length(setIndRun)
    
    idRun = setIndRun(iRun);
    matTS_movie(:, :, iRun) = tSeries_DFL(idRun).mnTS_C_raw(:, stimTiming_DFL(idRun).locCaFrame.locMovieOn:stimTiming_DFL(idRun).locCaFrame.locMovieOff); 
end


indNeuron = neuron.orderROIs('snr');
setCell = indNeuron;
figCheck = figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 880 315])
for iCell = 1:length(setCell) %1:size(matTS_movie, 1)
    idCell = setCell(iCell);
    figure(figCheck);
    SP(1) = subplot(2, 1, 1);
    imagesc(zscore(squeeze(matTS_movie(idCell, :, :)))')
    colormap(hot)
    ylabel(SP(1), 'Viewing')
    set(SP(1), 'CLim', [0 8]);
%     plot(squeeze(matTS_movie(iCell, :, :)))
    SP(2) = subplot(2, 1, 2);
    plot(zscore(squeeze(matTS_movie(idCell, :, :))))
    axis tight
    title(SP(1), sprintf('Cell #%d/%d (Cell ID: %d): movie %s', iCell, size(matTS_movie,1), tSeries_DFL(1).idNeuron_org(idCell), setMovie{iMovie}));

    set(SP(:), 'XTickLabel', 20:20:120)
    xlabel(SP(2), 'Time (s)');
    ylabel(SP(2), 'Norm. resp. (std)')
    
    F(iCell) = getframe(gcf);
    drawnow    
end

fname_movie = sprintf('%s_%s_DFL_%dCells_orderedSNR_tSeries_movie%s', dateSession, nameSubj, length(setCell), setMovie{iMovie});
writerObj = VideoWriter(fullfile(dirProcdata_session, fname_movie), 'MPEG-4');
writerObj.FrameRate = 1;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);



% %% selectivity calculation
% [sortedCond, indTrialCond] = sortrows(condMat);
% catMatAmp = cat(2, resultsTrial(indTrialCond).matAmp);
% catMnAmp = cat(2, resultsTrial(indTrialCond).mnAmp);
% figure
% plot(catMatAmp(5,:), 'o')
% figure
% plot(catMatAmp(14,:), 'o');  hold on
% line([110 209 319 429; 110 209 319 429], [-3 -3 -3 -3; 7 7 7 7], 'Color', 'm')
% 
% clear tempCat
% setCell = [5 14 37 38 46 47 48 54 76 77]; %1:size(tSeries_DFL_Cond, 1); %[5 14 37 38 46 47 48]; % [1 4 6];
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell); %4; %1;
%     for iCond = 1:5
%         tempCat{iCond} = cat(2, tSeries_DFL_Cond(idCell, :, iCond).matTS);
%     end
%     for iCond = 1:5
%         mnTS(:,iCond) = mean(tempCat{iCond}, 2);
%         steTS(:,iCond) = std(tempCat{iCond}, [], 2)./sqrt(size(tempCat{iCond}, 2));
%     end
%     
%     figure(100);
%     set(gcf, 'Position', [100 100 400 500])
%     plot([1:51].*0.1, mnTS, 'LineWidth', 2)
% %     legend('HF', 'MF', 'NO', 'FO', 'RD')
%     title(sprintf('Cell #%d', iCell))
%     axis tight
%     set(gca, 'Box', 'off', 'TickDIr', 'out', 'LineWidth', 2)
%     
%     
%     print(gcf, fullfile(dirFig, sprintf('%s_%s_ConditionMerge_Cell%d', dateSession, nameSubj, iCell)), '-depsc')
% 
% end
% 
% neuron.show_contours([], neuron.ids(indOrder_neuron(setCell)), neuron.PNR.*neuron.Cn, 'true');
% print(gcf, fullfile(dirFig, sprintf('%s_%s_cellContour_orderedSNR_Cell%s_label', dateSession, nameSubj, num2str(setCell))), '-depsc')
% 
% %% Bonus: show contours of cells
% addpath('/projects/parksh/_toolbox/CNMF_E/');
% cnmfe_setup;
% 
% nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
% dateSession = '20191113'; %'20191125';
% % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
% dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);
% d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D*'));
% load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); % man
% indOrder_neuron = neuron.orderROIs('snr');
% 
% neuron_b1=neuron.batches{1}.neuron;
% neuron_b2=neuron.batches{2}.neuron;
% 
% setCell = [5 14 37 38 46 47 48 54 76 77]; %
% neuron.ids(indOrder_neuron(setCell)) % in case of 20191113_Tabla's Sources2D results, when you order ROIs, the neuron.ids is different from neuron.batches.neuron
% 
% neuron_b1.show_contours([], [], neuron_b1.PNR.*neuron_b1.Cn, 'true');
% dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';
% print(gcf, fullfile(dirFig, '20191113_Tabla_cellContour_orderedSNR_77Cells'), '-depsc')
% 
% neuron_b1.show_contours([], neuron.ids(indOrder_neuron(setCell)), neuron_b1.PNR.*neuron_b1.Cn); %, 'true');
% print(gcf, fullfile(dirFig, sprintf('20191113_Tabla_cellContour_orderedSNR_Cell%s', num2str(setCell))), '-depsc')




% %% sort the timeseries
% [sortedCond, i] = sortrows(condMat);
% setCond = unique(sortedCond(:,1));
% 
% cond = cat(1, data.Condition(data.TrialError>1));
% [idCond, locTrial] = sort(cond);
% 
% setCond = unique(idCond);
% nCond = length(setCond);
% 
% % clear matTS_sorted %= [];
% matTS_sorted = [];
% for iCond = 1:nCond
%     curLoc = locTrial(idCond==setCond(iCond));
%     catCond = cat(3, matTS_trial{curLoc});
%     matTS_sorted = cat(4, matTS_sorted, catCond); % time x cell x trials x conditions
% end

% %% quick check across runs
% flagCheck = 1; %0; % if you want to go over cells to check quality across runs
% 
% if flagCheck
%     fig_check = figure; 
%     clear SP
%     
%     for iCell = 1:size(tSeries_DFL,1)
%         figure(fig_check); clf;
%         nPixel=10; % tSeries_DFL(iCell,1).nPixel;
%         
%         SP(1) = subplot(3,1,1);
%         %     plot(tSeries_DFL(iCell, 1).matTS');
%         
%         plot(tSeries_DFL(iCell,1).mnTS_C_raw, 'k-');
%         hold on
% %         line([1:length(tSeries_DFL(iCell,1).mnTS); 1:length(tSeries_DFL(iCell,1).mnTS)],...
% %             [tSeries_DFL(iCell,1).mnTS-tSeries_DFL(iCell,1).steTS; tSeries_DFL(iCell,1).mnTS+tSeries_DFL(iCell,1).steTS], 'Color', 'k')
%         title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 1, nPixel, max(tSeries_DFL(iCell,1).mnTS_C_raw)))
%         
%         SP(2) = subplot(3,1,2);
%         %     plot(tSeries_DFL(iCell, 2).matTS');
%         %     hold on
%         plot(tSeries_DFL(iCell,2).mnTS_C_raw, 'k-');
%         hold on
% %         line([1:length(tSeries_DFL(iCell,2).mnTS); 1:length(tSeries_DFL(iCell,2).mnTS)],...
% %             [tSeries_DFL(iCell,2).mnTS-tSeries_DFL(iCell,2).steTS; tSeries_DFL(iCell,2).mnTS+tSeries_DFL(iCell,2).steTS], 'Color', 'k')
%         title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 2, nPixel, max(tSeries_DFL(iCell,2).mnTS_C_raw)))
%         
%         SP(3) = subplot(3,1,3);
%         %     plot(tSeries_DFL(iCell, 3).matTS');
%         %     hold on
%         plot(tSeries_DFL(iCell,3).mnTS_C_raw, 'k-');
%         hold on
% %         line([1:length(tSeries_DFL(iCell,3).mnTS); 1:length(tSeries_DFL(iCell,3).mnTS)],...
% %             [tSeries_DFL(iCell,3).mnTS-tSeries_DFL(iCell,3).steTS; tSeries_DFL(iCell,3).mnTS+tSeries_DFL(iCell,3).steTS], 'Color', 'k')
%         title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 3, nPixel, max(tSeries_DFL(iCell,3).mnTS_C_raw)))
%         
%         axis(SP, 'tight')
%         input('')
%         
%     end
% end

