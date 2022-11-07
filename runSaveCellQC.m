% runSaveCellQC.m
%
% 2022/11/03 SHP
% - save the variance of noise (scale factor) and save the valid cell
% indices based on the SNR from the entire ts of the session
% 2022/10/25 SHP
% - add the procedure of selecting out spatially more localized cells at
% the level of 80% (infoCells(iSession).indCellValid_spatialCluster_0p2)
% - remove the part that exclude cells from Max's FOV ROI-based methods
%   : instead, just carrying that info here (infoCells(iSession).indCellValid_fov)
%   so that it is available with the spatial-localization-based valid cell indices
% 2022/06/01 SHP
% - load each session and save the SNR, PNR, number of cells, etc.
%

clear all; close all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/NIFVAULT/projects/parksh';
    dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
    dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
else % on virtual machine
    dirProjects = '/nifvault/projects/parksh';
    dirProcdata = '/nifvault/procdata/parksh';
    dirRawdata = '/nifvault/rawdata/parksh';
end

flagSaveFile = 1;


%% Get session info
setSubj = {'Tabla', 1; 'Max', 3};

for iSubj = 1:size(setSubj,1)
    
    nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
    FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
    [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    infoCells = struct([]);
    
    for iSession = 1:length(setDateSession)
        
        dateSession = setDateSession{iSession}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
        % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
        
        dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');
        
        
        %% Read source data
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        %% Cell quality check using SNR and PNR
        snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2); % signal variance divided by noise variance
        pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2); % peak amplitude divided by noise std
        varNoise = var(neuron.C_raw-neuron.C, 0, 2); % variance of noise
        
        [a, sortedID_snr] = sort(snrs, 'descend');
        [aa, sortedID_pnr] = sort(pnrs, 'descend');
        
        % % draw all the contours
        % neuron_b.show_contours([], [], imgFOV, 'true');
        
        imgFOV = neuron.Cn.*neuron.PNR;
        [center] = neuron.estCenter();
        
        thr = 0.2;
        Coor = neuron.get_contours(thr);
        [nrows, ncols] = cellfun(@size, Coor);
        locClustered = find(ncols > 2); % the ones with a reasonable localized contour at this threshold
        
        
        %     figure
        %     imagesc(imgFOV)
        %     hold on
        %     plot(center(indCell_highSNR,2), center(indCell_highSNR, 1), 'g.')
        %     text(center(indCell_highSNR,2)+1, center(indCell_highSNR,1), num2str([1:size(indCell_highSNR,1)]'), 'Color', 'g');
        %     axis off
        
        %% concatenate data into a struct
        infoCells(iSession).snrs = snrs;
        infoCells(iSession).pnrs = pnrs;
        infoCells(iSession).varNoise = varNoise;
        infoCells(iSession).snrs_mean = mean(snrs);
        infoCells(iSession).snr_crit = 0.3;
        infoCells(iSession).indCellValid_snr = find(snrs>infoCells(iSession).snr_crit);
%         infoCells(iSession).indCell_highSNR = find(snrs>mean(snrs));
%         infoCells(iSession).indCell_highPNR = find(snrs>mean(pnrs));
        infoCells(iSession).indCellValid_fov = (1:length(neuron.ids))';
        if strcmpi(nameSubj, 'max')
            load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
            infoCells(iSession).indCellValid_fov = indCell.validCell;
        end
        infoCells(iSession).coor_0p2 = Coor;
        infoCells(iSession).indCellValid_spatialCluster_0p2 = locClustered;
        %     infoCells(iSession).coor_0p5 = neuron.get_contours(0.5);
        
        
        % Basic info
        infoCells(iSession).dirSession = dirProcdata_session;
        infoCells(iSession).imgFOV = imgFOV;
        infoCells(iSession).cellCenter = center;
        infoCells(iSession).A = neuron.A;
        
        fprintf(1, '\n processed session #%d/%d...', iSession, length(setDateSession))
        
    end
    
    %% Save the data
    %% in main procdata FOV folder
    if flagSaveFile
        fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
            nameSubj, FOV_ID, nameSubj, FOV_ID));
        save(fname_cellQC, 'infoCells')
    end
    
end

%     %% try out DFL signal noise variance calculation
%     load(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'), 'tS_session');
%
%     % compute snr based only on the movie responses
%
%
%     for iCell = 1:size(tS_session(1).matTS_C_raw, 2)
%         aa = squeeze(tS_session(1).matTS_C_raw(:,iCell,:));
%         yyyRaw = aa(:);
%         aaC = squeeze(tS_session(1).matTS_C(:,iCell,:));
%         yyyC = aaC(:);
%         yyyN = yyyRaw-yyyC;
%
%         [pp, ff] = periodogram(yyyRaw, [], length(yyyRaw), 10);
%         locNoiseF = find(ff>2.5);
%         locSignalF = find(ff<1);
%         powerRaw_NoiseF(iCell,1) = sum(pp(locNoiseF));
%         powerRaw_SignalF(iCell,1) = sum(pp(locSignalF));
%
%         % compute other measures
%         snr_usingC(iCell,1) = var(yyyC)./var(yyyN);
%         pnr_usingC(iCell,1) = max(yyyC)./var(yyyN);
%         noise_usingC(iCell,1) = var(yyyN);
%
%     end


% %% playing
% for i = 1:length(infoCells)
%     indlow10 = [];
%     indlow10 = find(infoCells(i).pnrs<10);
%
%     [sortpnrs, indsort] = sort(infoCells(i).pnrs, 'descend');
%     indtop10 = indsort(1:10);
%
%     figure(100)
%     plot(infoCells(i).cellCenter(indtop10,2), infoCells(i).cellCenter(indtop10, 1), '.', 'MarkerSize', 15);
%     set(gca, 'YDir', 'reverse')
%     hold on;
%     input('')
% end

% %% plot movie tseries of potential same cell
% cellID = [22 16 14 17 19 18 nan nan nan 5 9 24]; % from Tabla FOV1
% nameSubj = 'Tabla'; %'Max'; %'Tabla';
% FOV_ID = 1; %3; %1;
% [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
%
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);
% tempMatTS1 = []; tempMatTS2 = [];
% for iS = 1:length(cellID)
%
%     if isnan(cellID(iS))
%         continue;
%     end
%
%     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
%     % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
%
%     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%
%     load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
%
%     figure(200);
%     subplot(2,1,1)
%     title('Mov 1')
%     plot(squeeze(tS_session(1).matTS_norm(:, cellID(iS), :)))
%     hold on
%
%     tempMatTS1 = cat(1, tempMatTS1, squeeze(tS_session(1).matTS_norm(:, cellID(iS), :))');
%
%     subplot(2,1,2)
%     title('Mov 2')
%     plot(squeeze(tS_session(2).matTS_norm(:, cellID(iS), :)))
%     hold on
%
%     tempMatTS2 = cat(1, tempMatTS2, squeeze(tS_session(2).matTS_norm(:, cellID(iS), :))');
% end


%     indCell_highSNR = find(snrs>mean(snrs));
%
%     for i = 1:124
%         figure(100);
%         subplot(2,1,1); cla;
%         plot(neuron.C_raw(sortedID_snr(i), 1000:2000))
%         hold on
%         plot(neuron.C(sortedID_snr(i), 1000:2000), 'r-')
%         axis tight
%         title(sprintf('Session %s: cell %d SNR = %2.2f', [nameSubj dateSession], sortedID_snr(i), a(i)))
%         subplot(2,1,2)
%         plot(neuron.S(sortedID_snr(i), 1000:2000))
%         axis tight
%         input('')
%     end
%
%
% load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
% load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
% load(fullfile(dirProcdata_session, 'RS_ts.mat'))
%     for i = 1:124
%         figure(200);
%         subplot(3,1,1); cla;
%         plot(tSeries_BPM(1).C_raw(sortedID_pnr(i), 1:1200))
%         hold on
%         plot(tSeries_BPM(1).C(sortedID_pnr(i), 1:1200), 'c-')
%         plot(tSeries_BPM(2).C_raw(sortedID_pnr(i), 1:1200), 'r-')
%         plot(tSeries_BPM(2).C(sortedID_pnr(i), 1:1200), 'm-')
%         axis tight
%         title(sprintf('%s BPM: cell %d SNR = %2.2f', [nameSubj dateSession], sortedID_pnr(i), a(i)))
%
%         subplot(3,1,2); cla;
%         plot(tSeries_DFL(4).C_raw(sortedID_pnr(i), :))
%         hold on
%         plot(tSeries_DFL(4).C(sortedID_pnr(i), :), 'c-')
%         plot(tSeries_DFL(5).C_raw(sortedID_pnr(i), :), 'r-')
%         hold on
%         plot(tSeries_DFL(5).C(sortedID_pnr(i), :), 'm-')
%         axis tight
%         title(sprintf('%s DFL: cell %d SNR = %2.2f', [nameSubj dateSession], sortedID_pnr(i), a(i)))
%
%         subplot(3,1,3); cla;
%         plot(tSeries_RS.C_raw(sortedID_pnr(i), 1001:2200))
%         hold on
%         plot(tSeries_RS.C(sortedID_pnr(i), 1001:2200), 'c-')
%         plot(tSeries_RS.C_raw(sortedID_pnr(i), 3001:4200), 'r-')
%         hold on
%         plot(tSeries_RS.C(sortedID_pnr(i), 3001:4200), 'm-')
%         axis tight
%         title(sprintf('%s RS: cell %d SNR = %2.2f', [nameSubj dateSession], sortedID_pnr(i), a(i)))
%
%         input('')
%     end




% get the contours and image field of view
% neuron_b = neuron.batches{1}.neuron;



%     % Generate cell location map within FOV
thr = 0.5; % the lower the smaller (more centralized) the contour
%     cellColor = [1 1 1];
%     widthContour = 1;
%     [d1,d2] = size(neuron.Cn);
%
%     figure;
%     imagesc(zeros(d1, d2)); % background
%     colormap(gray);
%     caxis([0 0.1]);
%     hold on;
%
%     CC = cell(size(neuron.A, 2),1);
%     CR = cell(size(neuron.A, 2),2);
%     for i = 1:size(neuron.A ,2)
%         A_temp = full(reshape(neuron.A(:,i),d1,d2));
%         A_temp = medfilt2(A_temp,[3,3]);
%         A_temp = A_temp(:);
%         [temp,ind] = sort(A_temp(:).^2,'ascend');
%         temp =  cumsum(temp);
%         ff = find(temp > (1-thr)*temp(end),1,'first');
%         if ~isempty(ff)
%             CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
%             fp = find(A_temp >= A_temp(ind(ff)));
%             [ii,jj] = ind2sub([d1,d2],fp);
%             CR{i,1} = [ii,jj]';
%             CR{i,2} = A_temp(fp)';
%         end
%         hold on;
%     end
%     axis off
%     title(sprintf('%s: %s', nameSubj, dateSession))
%     % %save
%     % print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');
%
%     % end
% %     %
%     Coor = neuron.get_contours(thr); % Coor = get_contours(obj, thr, ind_show); % ind_show: indices of cells you want to get contours
%     imgFOV = neuron.Cn.*neuron.PNR;
%
%     % % draw all the contours
%     % neuron_b.show_contours([], [], imgFOV, 'true');
%
%     figure
%     [center] = neuron.estCenter();
%     center = center(validIndCell, :);
%     imagesc(imgFOV)
%     hold on
%     plot(center(:,2), center(:, 1), 'r.')
%     text(center(:,2)+1, center(:,1), num2str([1:size(center,1)]'), 'Color', 'w');
%     axis off





% end



% % Using adjusted (ring background subtracted) time series from CNMFe
%     count = 0;
%     for iRun = 1:length(listRun_BPM)
%
%         % get the # of frames from each session without loading the matrix
%         matObj = matfile(fullfile(dirProcdata_session, d_file_imaging(iRun).name));
%         [d1 d2 nFrame] = size(matObj, 'Y');
%
%         for iCell = 1:size(neuron.A, 2)
% %             clear tempTS
% %             tempTS = matDFF(neuron.A(:,iCell)>0, :);
%
%             tSeries_BPM(iCell, iRun).idNeuron_org =neuron.ids(iCell);
% %             tSeries_BPM(iCell, iRun).matTS = tempTS;
% %             tSeries_BPM(iCell, iRun).nPixel = size(tempTS, 1);
%             tSeries_BPM(iCell, iRun).mnTS_C_raw = neuron.C_raw(iCell, count+1:count+nFrame);  %mean(tempTS);
%             tSeries_BPM(iCell, iRun).mnTS_C = neuron.C(iCell, count+1:count+nFrame);
% %             tSeries_BPM(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
%
%         end
%         count = count+nFrame;
%     end




% listRun_BPM = c{1}(contains(c{2}, 'BPM'));
% listRun_DFL = c{1}(contains(c{2}, 'DFL'));
% listRun_RS = c{1}(contains(c{2}, 'RS'));



% list of BPM imaging files
% listRun_BPM = {'103640', '104434', '105216'}'; % for 20191125_Tabla % eventually you want to retrieve this info from separate log file or something
% listRun_BPM = {'122312', '123102', '123609'}'; % for 20191125_Max
% listRun_BPM = {'123055', '123716', '124400', '125049', '125746', '130541'}'; % for 20191113_Tabla
% listRun_BPM = {'110653', '111411'}'; %{'105847', '110653', '111411'}'; %20191113_Max
