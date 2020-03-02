% temp_analyzeBPM.m
%
% quick & dirty analysis of BPM sessions
% just checking the data and try out bunch of stuff, before making the code
% neat and tidy

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
nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);
dirPreproc = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, dateSession, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

dirBHV = fullfile(dirRawdata, '/behavior/MonkeyLogic_Ca/'); %

% list of BPM imaging files
% listRun_BPM = {'103640', '104434', '105216'}'; % for 20191125_Tabla % eventually you want to retrieve this info from separate log file or something
% listRun_BPM = {'122312', '123102', '123609'}'; % for 20191125_Max
listRun_BPM = {'123055', '123716', '124400', '125049', '125746', '130541'}'; % for 20191113_Tabla
% listRun_BPM = {'110653', '111411'}'; %{'105847', '110653', '111411'}'; %20191113_Max

%% Copy DFF data from /_preproc  with renaming of "BPM_1" 
d_file_imaging = dir(fullfile(dirProcdata_session, 'BPM*_dFF.mat'));
if isempty(d_file_imaging)
    d_dFF = dir(fullfile(dirPreproc, '*_dFF.mat'));
    listDFF = {d_dFF.name}';
    for iRun = 1:length(listRun_BPM)
        indFile = contains(listDFF, listRun_BPM{iRun});
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(dirPreproc, listDFF{indFile}), fullfile(dirProcdata_session, sprintf('BPM_%d_%s_dFF.mat', iRun, listRun_BPM{iRun})))
    end
else
    fprintf(' *_dFF.mat files are already present in %s\n', dirProcdata_session)
end

%% Read source data
addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/')); 
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_BPM*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); % manually copied and renamed the file this time
% Deleted 5 neurons:  2, 48, 74, 75, 78, from Tabla's Sources2D file
% (Sources2D_02-Dec_10_50_16.mat)

indOrder_neuron = neuron.orderROIs('snr');
% matSource_org = neuron.reshape(neuron.A, 2); % 3D structure of spatial information of neurons
% nNeuron_org = size(matSource_org, 3);


%% Read each run's DFF and timing files
if ~exist(fullfile(dirProcdata_session, 'BPM_ts_tML.mat'), 'file')
    
    d_file_imaging = dir(fullfile(dirProcdata_session, 'BPM*_dFF.mat'));
%     d_file_timing = dir(fullfile(dirProcdata_session, 'BPM*_tML.mat'));
    
% % Getting tSeries from dFF files computed outside of CNMFe
%     for iRun = 1:length(listRun_BPM)
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
%             tSeries_BPM(iCell, iRun).idNeuron_org =neuron.ids(iCell);
%             tSeries_BPM(iCell, iRun).matTS = tempTS;
%             tSeries_BPM(iCell, iRun).nPixel = size(tempTS, 1);
%             tSeries_BPM(iCell, iRun).mnTS = mean(tempTS);
%             tSeries_BPM(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
%             
%         end        
%         clear Y_dFF matDFF
%         
%     end
    
    % Using adjusted (ring background subtracted) time series from CNMFe
    count = 0;
    for iRun = 1:length(listRun_BPM)
        
        % get the # of frames from each session without loading the matrix
        matObj = matfile(fullfile(dirProcdata_session, d_file_imaging(iRun).name));
        [d1 d2 nFrame] = size(matObj, 'Y');

        for iCell = 1:size(neuron.A, 2)
%             clear tempTS
%             tempTS = matDFF(neuron.A(:,iCell)>0, :);
            
            tSeries_BPM(iCell, iRun).idNeuron_org =neuron.ids(iCell);
%             tSeries_BPM(iCell, iRun).matTS = tempTS;
%             tSeries_BPM(iCell, iRun).nPixel = size(tempTS, 1);
            tSeries_BPM(iCell, iRun).mnTS_C_raw = neuron.C_raw(iCell, count+1:count+nFrame);  %mean(tempTS);
            tSeries_BPM(iCell, iRun).mnTS_C = neuron.C(iCell, count+1:count+nFrame);
%             tSeries_BPM(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
            
        end   
        count = count+nFrame;
    end
        
   resultsCNMFe.A = neuron.A;
   resultsCNMFe.options = neuron.options;
   resultsCNMFe.Cn = neuron.Cn;
   resultsCNMFe.PNR = neuron.PNR;
   resultsCNMFe.ids = neuron.ids;
   resultsCNMFe.indOrderNeuron = indOrder_neuron;
    
    d_file_timing = dir(fullfile(dirProcdata_session, 'BPM*_tML.mat'));
    
    for iRun = 1:length(listRun_BPM)
        
        load(fullfile(dirProcdata_session, d_file_timing(iRun).name), 'DataFile', 't_adj', 'stim')
        
        stimTiming_BPM(iRun).filename_org = DataFile;
        stimTiming_BPM(iRun).t_adj = t_adj;
        
        stimTiming_BPM(iRun).stim = stim;
        stimTiming_BPM(iRun).indValidTrial = find(stim.trialError>1);
        
        fs = 10; %eventually should be retrieved from xml file instead of hard-coding
        locTrialStart = floor(t_adj.trialStart./1000./(1/fs));
        locStimOn = floor(t_adj.stimOnset./1000./(1/fs));
        locBlankOn_afterStim = floor(t_adj.blankOnset_afterStim./1000./(1/fs));
        %     t_onset_adj = t_onset/1000-delay;
        %     locStimOn = floor(t_onset_adj./(1/fs));
        %     locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));
        
        stimTiming_BPM(iRun).locCaFrame.indTrial = find(stim.trialError>1);
        stimTiming_BPM(iRun).locCaFrame.locTrialStart = locTrialStart(stimTiming_BPM(iRun).locCaFrame.indTrial);
        stimTiming_BPM(iRun).locCaFrame.locStimOn = locStimOn;
        stimTiming_BPM(iRun).locCaFrame.locBlankOn_afterStim = locBlankOn_afterStim;
        
        
    end
    
    % interim save
    save(fullfile(dirProcdata_session, 'BPM_ts_tML.mat'), 'tSeries_BPM', 'stimTiming_BPM', 'resultsCNMFe')
    
else % if the file already exists
    load(fullfile(dirProcdata_session, 'BPM_ts_tML.mat'))
    fprintf(1, 'Time series data in %s is loaded to the workspace. \n', fullfile(dirProcdata_session, 'BPM_ts_tML.mat'))
end

%% quick check across runs
flagCheck = 1; %0; % if you want to go over cells to check quality across runs

if flagCheck
    fig_check = figure; 
    clear SP
    
    for iCell = 1:size(tSeries_BPM,1)
        figure(fig_check); clf;
        nPixel=10; % tSeries_BPM(iCell,1).nPixel;
        
        SP(1) = subplot(3,1,1);
        %     plot(tSeries_BPM(iCell, 1).matTS');
        
        plot(tSeries_BPM(iCell,1).mnTS_C_raw, 'k-');
        hold on
%         line([1:length(tSeries_BPM(iCell,1).mnTS); 1:length(tSeries_BPM(iCell,1).mnTS)],...
%             [tSeries_BPM(iCell,1).mnTS-tSeries_BPM(iCell,1).steTS; tSeries_BPM(iCell,1).mnTS+tSeries_BPM(iCell,1).steTS], 'Color', 'k')
        title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 1, nPixel, max(tSeries_BPM(iCell,1).mnTS_C_raw)))
        
        SP(2) = subplot(3,1,2);
        %     plot(tSeries_BPM(iCell, 2).matTS');
        %     hold on
        plot(tSeries_BPM(iCell,2).mnTS_C_raw, 'k-');
        hold on
%         line([1:length(tSeries_BPM(iCell,2).mnTS); 1:length(tSeries_BPM(iCell,2).mnTS)],...
%             [tSeries_BPM(iCell,2).mnTS-tSeries_BPM(iCell,2).steTS; tSeries_BPM(iCell,2).mnTS+tSeries_BPM(iCell,2).steTS], 'Color', 'k')
        title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 2, nPixel, max(tSeries_BPM(iCell,2).mnTS_C_raw)))
        
        SP(3) = subplot(3,1,3);
        %     plot(tSeries_BPM(iCell, 3).matTS');
        %     hold on
        plot(tSeries_BPM(iCell,3).mnTS_C_raw, 'k-');
        hold on
%         line([1:length(tSeries_BPM(iCell,3).mnTS); 1:length(tSeries_BPM(iCell,3).mnTS)],...
%             [tSeries_BPM(iCell,3).mnTS-tSeries_BPM(iCell,3).steTS; tSeries_BPM(iCell,3).mnTS+tSeries_BPM(iCell,3).steTS], 'Color', 'k')
        title(sprintf('Cell #%d, Run #%d: nPixel: %d, maxDFF: %2.2f', iCell, 3, nPixel, max(tSeries_BPM(iCell,3).mnTS_C_raw)))
        
        axis(SP, 'tight')
        input('')
        
    end
end





%% sort the timeseries for each cell and each condition
for iRun = 1:size(tSeries_BPM, 2)
    
%     indTrial = stimTiming_BPM(iRun).locCaFrame.indTrial;
    indValidTrial = find(stimTiming_BPM(iRun).locCaFrame.locTrialStart>0);
    condMat = stimTiming_BPM(iRun).stim.condMat(stimTiming_BPM(iRun).locCaFrame.indTrial(indValidTrial),:);
    
    [sortedCond, indTrialCond] = sortrows(condMat);
    setCond = unique(sortedCond(:,1));
    
    matTS = cat(1, tSeries_BPM(:,iRun).mnTS_C_raw)'; %mnTS_C_raw)'; %cat(1, tSeries_BPM(:,iRun).mnTS)';
    
    clear matTS_trial
    for iTrial = 1:length(indValidTrial)
        matTS_trial(:,:,iTrial) = matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))-10:stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+40, :); %locEnd(iTrial), :); %locOnset_end(iTrial), :);
    end
    
    clear resultsTrial
    for iTrial = 1:length(indValidTrial)
        matAmp = [];
        win_ms = [500 1500]; % time window of amplitude calculation (in ms: zero is stim Onset) 
        fs = 10;
        win_frame = win_ms./(1000/fs);
        matAmp = matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(1): stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+win_frame(2), :)';
                
        resultsTrial(iTrial).window_ms = win_ms;
        resultsTrial(iTrial).win_frame = win_frame;
        resultsTrial(iTrial).matAmp = matAmp;
        resultsTrial(iTrial).mnAmp = mean(matAmp, 2);
    end
        
    resultsBPM(iRun).resultsTrial = resultsTrial;
        %matTS(stimTiming_BPM(iRun).locCaFrame.locStimOn(indValidTrial(iTrial))+50
    % matTS_trial = cat(3, sortTrial.matTS);
    % clear sortTrial
    
    for iCell = 1:size(matTS_trial, 2)
        
        for iCond = 1:length(setCond)
            tSeries_BPM_Cond(iCell, iRun, iCond).condMat = condMat(indTrialCond(sortedCond(:,1)==setCond(iCond)),:);
            tSeries_BPM_Cond(iCell, iRun, iCond).indTrial = indTrialCond(sortedCond(:,1)==setCond(iCond));
            tSeries_BPM_Cond(iCell, iRun, iCond).matTS = squeeze(matTS_trial(:, iCell, tSeries_BPM_Cond(iCell, iRun, iCond).indTrial)); % matTS_trial(:,:,sortCond(iCond).indTrial);
        end
        
    end
    
end


fig_selectivity = figure;

for iCell = 1:size(tSeries_BPM_Cond, 1)
    fprintf(1, 'Cell #%d: ', iCell)
    
%     iRun = 1;
    for iRun = 1:size(tSeries_BPM_Cond,2) 
%         curTS_catCond = [];
%         curTS_catCond = cat(2, tSeries_BPM_Cond(iCell, iRun, :).matTS);
        
        for iCond = 1:5
            figure(fig_selectivity)
            SP(5*(iRun-1)+iCond) = subplot(size(tSeries_BPM_Cond,2), 5, 5*(iRun-1)+iCond);
            plot(tSeries_BPM_Cond(iCell, iRun, iCond).matTS)
%             plot(curTS_catCond(:,:,iCond))
            ylim([-3 15])
        end
        
        title(SP(3), sprintf('Cell #%d', iCell))
        
    end
    
    input('')
end

%% selectivity calculation
[sortedCond, indTrialCond] = sortrows(condMat);
catMatAmp = cat(2, resultsTrial(indTrialCond).matAmp);
catMnAmp = cat(2, resultsTrial(indTrialCond).mnAmp);
figure
plot(catMatAmp(5,:), 'o')
figure
plot(catMatAmp(14,:), 'o');  hold on
line([110 209 319 429; 110 209 319 429], [-3 -3 -3 -3; 7 7 7 7], 'Color', 'm')

clear tempCat
setCell = [5 14 37 38 46 47 48 54 76 77]; %1:size(tSeries_BPM_Cond, 1); %[5 14 37 38 46 47 48]; % [1 4 6];
for iCell = 1:length(setCell)
    idCell = setCell(iCell); %4; %1;
    for iCond = 1:5
        tempCat{iCond} = cat(2, tSeries_BPM_Cond(idCell, :, iCond).matTS);
    end
    for iCond = 1:5
        mnTS(:,iCond) = mean(tempCat{iCond}, 2);
        steTS(:,iCond) = std(tempCat{iCond}, [], 2)./sqrt(size(tempCat{iCond}, 2));
    end
    
    figure(100);
    set(gcf, 'Position', [100 100 400 500])
    plot([1:51].*0.1, mnTS, 'LineWidth', 2)
%     legend('HF', 'MF', 'NO', 'FO', 'RD')
    title(sprintf('Cell #%d', iCell))
    axis tight
    set(gca, 'Box', 'off', 'TickDIr', 'out', 'LineWidth', 2)
    
%     input('')
    print(gcf, fullfile(dirFig, sprintf('%s_%s_ConditionMerge_Cell%d', dateSession, nameSubj, iCell)), '-depsc')

end

neuron.show_contours([], neuron.ids(indOrder_neuron(setCell)), neuron.PNR.*neuron.Cn, 'true');
print(gcf, fullfile(dirFig, sprintf('%s_%s_cellContour_orderedSNR_Cell%s_label', dateSession, nameSubj, num2str(setCell))), '-depsc')

%% Bonus: show contours of cells
addpath('/projects/parksh/_toolbox/CNMF_E/');
cnmfe_setup;

nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D*'));
load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); % man
indOrder_neuron = neuron.orderROIs('snr');

neuron_b1=neuron.batches{1}.neuron;
neuron_b2=neuron.batches{2}.neuron;

setCell = [5 14 37 38 46 47 48 54 76 77]; %
neuron.ids(indOrder_neuron(setCell)) % in case of 20191113_Tabla's Sources2D results, when you order ROIs, the neuron.ids is different from neuron.batches.neuron

neuron_b1.show_contours([], [], neuron_b1.PNR.*neuron_b1.Cn, 'true');
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';
print(gcf, fullfile(dirFig, '20191113_Tabla_cellContour_orderedSNR_77Cells'), '-depsc')

neuron_b1.show_contours([], neuron.ids(indOrder_neuron(setCell)), neuron_b1.PNR.*neuron_b1.Cn); %, 'true');
print(gcf, fullfile(dirFig, sprintf('20191113_Tabla_cellContour_orderedSNR_Cell%s', num2str(setCell))), '-depsc')




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

