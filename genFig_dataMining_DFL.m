% genFig_dataMining_DFL.m
%
% 2022/11/01 SHP: started new by applying longitudinally registered cells,
% instead of looking at session to session data
% 2021/12/21 SHP: working on the part that selects good cells
% 2021/09/08 SHP
% Digging DFL data to find out something
% started from "genFig_functionalMapFOV.m"

clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/PROJECTS/parksh';
        dirProcdata = '/Volumes/NIFVAULT/PROCDATA/parksh';
        dirRawdata = '/Volumes/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');


%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);


%% load saved files
% cells pooled across days
fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_stack, 'cellIDAcrossDay'); %, 'stackCellCenter')

% cell quality info 
fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_cellQC, 'infoCells')

% translational shift across days
fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));  
load(fname_shifts, 'shifts')

% aligned cells TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')


%%
indCellValid = find(cat(1, cellTS.nTrial1_total)>8); % cells that have more than 8 trials for movie 1

for iCell = 1:length(indCellValid)   
    
    matAvgTS1(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie1)'; % do I want to use this z-scored ts?
    matAvgTS2(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie2)'; %

end



%%
% clear all;
% 
% ss = pwd;
% if ~isempty(strfind(ss, 'Volume')) % if it's local
%     dirProjects = '/Volumes/NIFVAULT/projects/parksh';
%     dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
%     dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
% else % on virtual machine
%     dirProjects = '/nifvault/projects/parksh';
%     dirProcdata = '/nifvault/procdata/parksh';
%     dirRawdata = '/nifvault/rawdata/parksh';
% end
% 
% % setNameSubj = {'Tabla', 'Max'};
% flagSavePPTX = 0; %1;
% 
% % get session info
% nameSubj = 'Tabla';
% FOV_ID = 1;
% [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);
% 
% for iSession = 1:nSession
% % iSession = 1; 
% dateSession = setDateSession{iSession}; %'20191113'; %setDateSession{iSession};
% 
% dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% dirPreproc = fullfile(dirProcdata_session, '_preproc');
% 
% dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
% 
% 
% %% Read source data
% addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
% cnmfe_setup;
% d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
% 
% load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
% 
% validIndCell = [];
% validIndCell(:,1) = 1:length(neuron.ids);
% if strcmpi(nameSubj, 'max')
%     load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
%     validIndCell = indCell.validCell;
% end
% 
% 
% % get the contours and image field of view
% % neuron_b = neuron.batches{1}.neuron;
% 
% % Generate cell location map within FOV
% thr = 0.5; % the lower the smaller (more centralized) the contour
% cellColor = [1 1 1];
% widthContour = 1;
% [d1,d2] = size(neuron.Cn);
% 
% figure;
% imagesc(zeros(d1, d2)); % background
% colormap(gray);
% caxis([0 0.1]);
% hold on;
% 
% CC = cell(size(neuron.A, 2),1);
% CR = cell(size(neuron.A, 2),2);
% for i = 1:size(neuron.A ,2)
%     A_temp = full(reshape(neuron.A(:,i),d1,d2));
%     A_temp = medfilt2(A_temp,[3,3]);
%     A_temp = A_temp(:);
%     [temp,ind] = sort(A_temp(:).^2,'ascend');
%     temp =  cumsum(temp);
%     ff = find(temp > (1-thr)*temp(end),1,'first');
%     if ~isempty(ff)
%         CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
%         fp = find(A_temp >= A_temp(ind(ff)));
%         [ii,jj] = ind2sub([d1,d2],fp);
%         CR{i,1} = [ii,jj]';
%         CR{i,2} = A_temp(fp)';
%     end
%     hold on;
% end
% axis off
% title(sprintf('%s: %s', nameSubj, dateSession))
% % %save
% % print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');
% 
% % end
% %
% Coor = neuron.get_contours(thr); % Coor = get_contours(obj, thr, ind_show); % ind_show: indices of cells you want to get contours
% imgFOV = neuron.Cn.*neuron.PNR;
% 
% % % draw all the contours
% % neuron_b.show_contours([], [], imgFOV, 'true');
% 
% figure
% [center] = neuron.estCenter();
% center = center(validIndCell, :);
% imagesc(imgFOV)
% hold on
% plot(center(:,2), center(:, 1), 'r.')
% text(center(:,2)+1, center(:,1), num2str([1:size(center,1)]'), 'Color', 'w');
% axis off
% 
% 
% %% Cell quality check using SNR and PNR
% snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);
% pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2);
% 
% snrs = snrs(validIndCell);
% pnrs = pnrs(validIndCell);
% 
% [a, sortedID_snr] = sort(snrs, 'descend');
% [aa, sortedID_pnr] = sort(pnrs, 'descend');
% 
% % SNR distribution 
% figure
% histogram(snrs, 30)
% median(snrs)
% set(0, 'defaultfigurecolor', [1 1 1])
% xlabel('SNR')
% title(sprintf('%s: %s', nameSubj, dateSession))
% line([median(snrs) median(snrs)], get(gca, 'ylim'), 'Color', 'r')
% 
% % Location of a particular set cells
% setCells_rank = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr); %1:10
% figure;
% imagesc(imgFOV);
% colormap(gray);
% hold on
% plot(center(sortedID_snr(setCells_rank), 2), center(sortedID_snr(setCells_rank), 1), 'r.')
% hold on
% text(center(sortedID_snr(setCells_rank), 2)+1, center(sortedID_snr(setCells_rank), 1)+1, num2str([setCells_rank]'), 'Color', 'r');
% text(center(sortedID_pnr(setCells_rank), 2)+3, center(sortedID_pnr(setCells_rank), 1)+3, num2str([setCells_rank]'), 'Color', 'g');
% plot(center(sortedID_pnr(setCells_rank), 2), center(sortedID_pnr(setCells_rank), 1), 'g.')
% axis off
% 
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(100);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_snr(iCell)), 1:nTime)+10*(iCell-1), '-')
% tY(iCell) = mean(neuron.C_raw(validIndCell(sortedID_snr(iCell)), 1:nTime)+10*(iCell-1));
% end
% title('SNR-based sorting')
% set(gca, 'YTick', tY, 'YTickLabel', 1:nCell, 'TickDir', 'out')
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% xlabel('Time (s)')
% ylabel('Cell order')
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(200);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+10*(iCell-1), '-')
% tY(iCell) = mean(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+10*(iCell-1));
% end
% title('PNR-based sorting')
% set(gca, 'YTick', tY, 'YTickLabel', 1:nCell, 'TickDir', 'out')
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% xlabel('Time (s)')
% ylabel('Cell order')
% 
% 
% tY=[];
% nCell = 10;
% nTime = 1000;
% for iCell =1:nCell
% figure(200);
% hold on;
% plot(neuron.C_raw(validIndCell(sortedID_pnr(iCell)), 1:nTime)+3*(iCell-1), '-', 'LineWidth', 2)
% end
% set(gca, 'XTick', 200:200:nTime,  'XTickLabel', [200:200:nTime]./10);
% 
% 
% 
% %% Movie-driven signal
% load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
% 
% % snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);
% % pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2);
% % 
% % [a, sortedID_snr] = sort(snrs, 'descend');
% % [aa, sortedID_pnr] = sort(pnrs, 'descend');
% 
% setCell = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr); %1:10; %11:20; %10; 
% tY = []; ttY = [];
% figMovie_snr = figure;
% set(figMovie_snr, 'Position', [675    31   660   930], 'name', sprintf('%s %s: PNR based sorting', nameSubj, dateSession))
% SP(1) = subplot(1,2,1);
% SP(2) = subplot(1,2,2);
% hold(SP(:), 'on');
% 
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell);
%     figure(figMovie_snr);
%     plot(SP(1), squeeze(tS_session(1).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1), '-'); 
%     tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1)));
%     
%     plot(SP(2), squeeze(tS_session(2).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1), '-'); 
%     ttY(iCell) = mean(mean(squeeze(tS_session(2).matTS(:, validIndCell(sortedID_snr(idCell)), :))+10*(iCell-1)));
% end
% set(SP, 'YTick', tY, 'YTickLabel', setCell, 'TickDir', 'out')
% set(SP, 'XTick', 200:200:1200,  'XTickLabel', [200:200:1200]./10);
% axis(SP, 'tight')
% xlabel(SP, 'Time (s)')
% ylabel(SP(1), 'Cell')
% title(SP(1), 'Movie 1')
% title(SP(2), 'Movie 2')
% 
% %pnr
% setCell = 1:10; %length(sortedID_pnr)-9:length(sortedID_pnr);  %1:10; %11:20; %10; 
% tY = []; ttY = [];
% figMovie_pnr = figure;
% set(figMovie_pnr, 'Position', [675    31   660   930], 'name', sprintf('%s %s: PNR based sorting', nameSubj, dateSession))
% SP(1) = subplot(1,2,1);
% SP(2) = subplot(1,2,2);
% hold(SP(:), 'on');
% 
% scalefac = 10;
% for iCell = 1:length(setCell)
%     idCell = setCell(iCell);
%     figure(figMovie_pnr);
%     plot(SP(1), squeeze(tS_session(1).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1), '-'); 
%     tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1)));
%     
%     plot(SP(2), squeeze(tS_session(2).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1), '-'); 
%     ttY(iCell) = mean(mean(squeeze(tS_session(2).matTS(:, validIndCell(sortedID_pnr(idCell)), :))+scalefac*(iCell-1)));
% end
% set(SP, 'YTick', tY, 'YTickLabel', setCell, 'TickDir', 'out')
% set(SP, 'XTick', 200:200:1200,  'XTickLabel', [200:200:1200]./10);
% axis(SP, 'tight')
% xlabel(SP, 'Time (s)')
% ylabel(SP(1), 'Cell')
% title(SP(1), 'Movie 1')
% title(SP(2), 'Movie 2')
% 
% 
% %% pupil size change
% % part of eye data is not great. decided to focus on the later half of the
% % first movie, which contains body motion, object motion, face
% 
%     % get session info
%     [infoSession, opts] = readInfoSession(nameSubj);
%     S = table2struct(infoSession);
%     
%     % setExpName = {S.ExpName}';
%     setMLFilename = {S.MLFilename}';
%     
%     indDFLRuns = contains(setMLFilename, 'DFL') & cat(1, S.flagPreproc) > 0 & contains({S.stimulus}', 'set1_1'); %% containing "DFL" in filename AND flagPreproc value of 1
%     setFilename = setMLFilename(indDFLRuns);
%     
%     % Tabla
%     % 191113: eye signal not good (another range of x/y/pupil present during the first half)
%     % 191114 & 191118: eye signal good. Length of data points are 6-7ms shorter than 2min
%     % except the last run of 191118
%     % 191119: starting to lose x signal + y & pupil is bad in general
%     % 191120: x is lost, y& pupil was good in the first run, but not in the
%     % others
%     % 191121: x is lost, length of data points are 6-7ms shorter than 2min, y
%     % & pupil bad in some runs
%    % 191125: x is lost, for the first file signals look okay & length is fine 
%    
%     
%     for iFile = 1:length(setFilename)
%         filename = strcat(setFilename{iFile}, '.bhv2');
%         
%         dateSession = filename(1:6);
%         
%         if str2num(dateSession) < 191121
%             dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
%         else
%             dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
%         end
%         
%         % filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have
%         
%         %% Read the file
%         data = mlread(fullfile(dirBHV, filename)); % mlread(filename);
%         
%         %% eye data during stimulus on 
%         % Event Code Numbers & Names : TASK_START = 10; FP_ON = 20;
%         % WAIT_FOR_TR = 30; MOVIE_ON = 40; REWARD = 90; 
%         % TRIG onset = 900; TRIG offset = 990; (TTL from ML to Inscopix DAQ On & Off)
%         locStimOn = find(data.BehavioralCodes.CodeNumbers == 40);
%         time_stimOn = floor(data.BehavioralCodes.CodeTimes(locStimOn))
% %         data.AnalogData
%         
%         figure;
%         plot(data.AnalogData.Eye, '.')
%         line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
%         title(sprintf('%s: x & y gaze', filename))
%         axis tight
%         xlabel('Time')
%         ylabel('Voltage')
%         legend('X', 'Y')
%         
%         figure;
%         plot(data.AnalogData.General.Gen1, 'g.')
%         line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
%         title(sprintf('%s: pupil size change', filename))
%         axis tight
%         xlabel('Time')
%         ylabel('Voltage')
%         
%         input('')
%     end
% tempP = data.AnalogData.General.Gen1((time_stimOn:end));
% figure
% plot(tempP)
% plot(tempP(60001:120000))
% for iCell = 1:length(setCell)
% figure(200); hold on;
% idCell = setCell(iCell);
% plot(squeeze(tS_session(1).matTS(:, sortedID_pnr(idCell), 1))+scalefac*(iCell-1), '-');
% tY(iCell) = mean(mean(squeeze(tS_session(1).matTS(:, sortedID_pnr(idCell), 1))+scalefac*(iCell-1)));
% end
% axis tight
% xlim([600 1200])


%% OLD CODE: 1) consistency across trials, 2) clustering of averaged movie-driven responses
% ts = struct([]);
% iMovie = 1;
% for iCell = 1:size(tS_session(iMovie).avgTS, 2)
%     
%     curMatTS = squeeze(tS_session(iMovie).matTS_norm(:,iCell,:));
%     avgMatTS = tS_session(iMovie).avgTS_norm(:,iCell);
%     steMatTS = std(curMatTS, [], 2)./sqrt(size(curMatTS, 2)-1);
%     
%     ts(iCell).avgMatTS = avgMatTS;
%     ts(iCell).steMatTS = steMatTS;
%     
%     % figure(100);clf;
%     % subplot(2,1,1)
%     % plot(curMatTS); axis tight
%     % title(sprintf('Cell #%d/%d', iCell, size(tS_session(iMovie).avgTS, 2)))
%     % subplot(2,1,2)
%     % plot(avgMatTS, 'k'); axis tight
%     % line(repmat(1:length(avgMatTS), 2, 1), cat(2, avgMatTS+steMatTS, avgMatTS-steMatTS)', 'Color', 'c')
%     
%     % input('')
% end
% 
% % %% Consistency across trials
% catAvgMatTS = cat(2, ts.avgMatTS); % driven activity, averaged across trials
% catSteMatTS = cat(2, ts.steMatTS);
% 
% [coeff, score, latent, tsquared, explained] = pca(catAvgMatTS');
% 
% k = 4;
% [IDXdfl, C, SUMD] = kmeans(catAvgMatTS', k, 'Distance', 'correlation');
% [sortedIDXdfl, indCelldfl] = sort(IDXdfl);
% clear indCell_sort
% for iType = 1:k
%     indCell_sort{iType} = indCelldfl(sortedIDXdfl==iType);
% end
% 
% % Plotting
% fig_summary_DFL = figure;
% set(gcf,  'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1085 750])
% clear sp
% 
% % 1. averaged amplitude for each condition
% figure(fig_summary_DFL);
% sp(1) = subplot('Position', [0.2 0.65 0.75 0.3]);
% imagesc(catAvgMatTS(:, indCelldfl)');
% colormap(hot);
% set(sp(1), 'CLim', [-1 1].*4)
% set(sp(1), 'YTick', find(diff(sortedIDXdfl)>0), 'YTickLabel', [])
% set(sp(1), 'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
% set(sp(1), 'TickDir', 'out')
% box off
% colorbar;
% title(sprintf('%s %s: averaged response to movie %s', nameSubj, dateSession, tS_session(iMovie).idStim))
% ylabel('Cells (sorted)')
% xlabel('Time (s)')
% 
% % 2. Clustering results on 2-d PC space
% figure(fig_summary_DFL);
% sp(2) = subplot('Position', [0.1 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% for iType = 1:k
%         plot(score(indCell_sort{iType}, 1), score(indCell_sort{iType}, 2), 'o', 'MarkerFaceColor', cMap_sort(iType, :));
%         hold on;
% end
% legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Location', 'best')
% xlabel(sprintf('PC 1: explained %2.2f %% var', explained(1)))
% ylabel(sprintf('PC 2: explained %2.2f %% var', explained(2)))
% set(gca, 'TickDir', 'out')
% box off
% axis square
% title('Clustering based on movie response on PC space')
% 
% % 3. Clustering results on imaging field of view
% figure(fig_summary_DFL);
% sp(3) = subplot('Position', [0.55 0.1 0.4 0.4]);
% cMap_sort = hsv(k);
% 
% imagesc(imgFOV); 
% colormap(sp(3), gray);
% hold on;
% for iType = 1:k
%     for iC = 1:size(indCell_sort{iType}, 1)
%         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
%     end
% end
% axis off
% 
% % % for each run
% % fig_map = figure;
% % for iRun = 1:length(tSeries_DFL)
% % matTSnorm = zscore(tSeries_DFL(iRun).C_raw');
% % 
% % k = 5;
% % [IDX, C, SUMD] = kmeans(matTSnorm', k, 'Distance', 'correlation');
% % [sortedIDX, indCell] = sort(IDX);
% % 
% % 
% % clear indCell_sort
% % for iType = 1:k
% %     indCell_sort{iType} = indCell(sortedIDX==iType);
% % end
% % 
% % cMap_sort = hsv(k);
% % figure(fig_map);
% % subplot(1,length(tSeries_DFL),iRun);
% % imagesc(imgFOV); colormap(gray);
% % axis off
% % for iType = 1:k
% %     for iC = 1:size(indCell_sort{iType}, 1)
% %         figure(fig_map);
% %         hold on;
% %         plot(Coor{indCell_sort{iType}(iC, 1)}(1,:), Coor{indCell_sort{iType}(iC, 1)}(2,:), '.', 'Color', cMap_sort(iType, :));
% %     end
% % end
% % title(sprintf('DFL Run #%d', iRun))
% % 
% % end
% 
% % end
% 
% % if flagSavePPTX
% %     % save figures
% %     addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
% %     addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
% %     
% %     fname_pptx = sprintf('%s_ClusteringBPMDFL', nameSubj); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
% %     exportFigsToPPTX(fname_pptx);
% %     
% % %     switch lower(nameSubj)
% % %         case 'tabla'
% % %             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
% % %         case 'max'
% % %             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
% % %     end
% %     movefile('./*.pptx', dirFig);
% % end
% 
% % end
% 
% 
% 


