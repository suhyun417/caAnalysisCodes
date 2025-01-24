% genFig_exampleFOV_TS.m
% 
% 2024/05/23 SHP
% Generate figure of example FOV and calcium traces 


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

iSubj = 1; %2; %1; %2; %1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%%
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

% % aligned cells TS and spatial info
% fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
% load(fname_caTSFOV_DFL, 'cellTS', 'cellPix')
% cellTS_DFL = cellTS;
% clear cellTS
% 
% fname_caTSFOV_BPM = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
% load(fname_caTSFOV_BPM, 'cellTS')
% cellTS_BPM = cellTS;
% clear cellTS


%% FOV
iSession = 1; %:length(setDateSession) %;
dateSession = setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

% dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

load(sprintf('%s/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', directory.dirProcdata, nameSubj, dateSession))
% load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

% Read source data
addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));

load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));

% get the contours and image field of view
% thr = 0.3; % the lower the smaller (more centralized) the contour
% Coor = neuron.get_contours(thr);
% imgFOV = neuron.Cn.*neuron.PNR;
% 
% % draw all the contours
% neuron.show_contours([], [], imgFOV, 'true');

% draw contours for each session
% figure;
% subplot('Position', [0 0 1 1]);
% imagesc(neuron.Cn.*neuron.PNR); colormap(gray);
% hold on;

figure;
subplot('Position', [0 0 1 1]);
image(ones(size(neuron.Cn)).*255);
colormap(gray);
hold on

thr = 0.5; %0.2;
cellColor = ones(1,3).*0.75; %'k'; %cMap_line(iSession, :); %'m'; %'c'; 
widthContour = 2;
[d1,d2] = size(neuron.Cn);
indCellValid_session = cellIDAcrossDay(~isnan(cellIDAcrossDay(:,iSession)), iSession);

CC = cell(length(indCellValid_session),1);
CR = cell(length(indCellValid_session),2);
% cmap_cell = cool(size(neuron.A, 2));
for iCC = 1:length(indCellValid_session)
    i = indCellValid_session(iCC); %sortedIndCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
truesize;
set(gcf, 'PaperPositionMode', 'auto')

% add scale bar of 100um (=32px)
hold on;
line([1 33], [265 265], 'Color', 'k', 'LineWidth', 5)

% print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_gray_scaleBar_r0', ...
%     nameSubj, FOV_ID, thr*10, iSession)), '-r0', '-depsc');


%% Example cell TS: select various examples depending on BPM>DFL, BPM<DFL, BPM=DFL 
% DFL_ts.mat,BPM_ts.mat, RS_ts.mat: raw Ca TS (not-normalized, just cut for each run from concatenated source extraction)
% BPM>DFL: Tabla session 1 cell 40, 9, 73, 82
% BPM<DFL: Tabla session 1 cell 21, 4, 35, 41 
% BPM=DFL: Tabla session 1 cell 79, 36, 23, 29
% select the cells 

% first check the PNRs consistency across runs
load(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'))
pnrs_dfl_set = squeeze(max(tS_session(1).matTS_C)./std(tS_session(1).matTS_C_raw-tS_session(1).matTS_C)); % for all the runs

load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
for iRun = 1:length(tSeries_BPM)
pnrs_bpm_set(:,iRun) = max(tSeries_BPM(1).C, [], 2)./std(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2); % p
end

load(fullfile(dirProcdata_session, 'BPM_ts_tML.mat'), 'tS_run')

%%% 1. draw all cells together, focus on MOV
% set of the example cells 
setCell = [9 73 82 29 36 79 41 4 35 21]; % changed order to show the larger movie response at the top

% prepare the cell traces
matTS_BPM = [];
for iCell = 1:length(setCell)
    matTS_BPM(:, iCell) = cat(1, tS_run(1).tS_trial(setCell(iCell),:).matTS); % this is C_raw
end

iRun = 1; iMovie = 1;
matTS_DFL = [];
matTS_DFL = tS_session(iMovie).matTS_C_raw(:, setCell, iRun);

load(fullfile(dirProcdata_session, 'RS_ts.mat'))
matTS_RS = [];
matTS_RS = tSeries_RS(1).C_raw(setCell, :)';


cMap_lineColor = [102 194 165; 252 141 98; 141 106 203]./255;
cMap_lineColor2 = [27 158 119; 217 95 2; 117 112 179]./255;

% draw all the TS together: this time MOV-IMG-RS, and cell with larger MOV
% responses come first
widthBreak = 50;
matTS_grand = [];
matTS_grand = cat(1, matTS_DFL(1:1200,:), NaN(widthBreak, length(setCell)),...
    matTS_BPM(1:1200, :), NaN(widthBreak, length(setCell)),...
    matTS_RS(2100:2100+1199, :));

fig_grand = figure;
ylim = [2 36];
BPM_stimOn = [11:46:1200]+1200+widthBreak;
line(repmat(BPM_stimOn, 2, 1), repmat(ylim', 1, length(BPM_stimOn)), 'Color', ones(1,3).*0.5)
hold on;
DFL_stimOn = [200:200:1200];
line(repmat(DFL_stimOn, 2, 1), repmat(ylim', 1, length(DFL_stimOn)), 'Color', ones(1,3).*0.5)

p = plot(matTS_grand+repmat([1:length(setCell)].*3, size(matTS_grand, 1), 1), 'LineWidth', 2);
set(p(1:3), 'Color', cMap_lineColor(1,:))
set(p(4:6), 'Color', cMap_lineColor(2,:))
set(p(7:10), 'Color', cMap_lineColor(3,:))
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'k', 'LineWidth', 2)
axis tight
set(gca, 'YLim', [2 36])

% BPM_stimOn = [11:46:1200]+1200+widthBreak;
% line(repmat(BPM_stimOn, 2, 1), repmat(get(gca, 'ylim')', 1, length(BPM_stimOn)), 'Color', ones(1,3).*0.5)
% DFL_stimOn = [200:200:1200];
% line(repmat(DFL_stimOn, 2, 1), repmat(get(gca, 'ylim')', 1, length(DFL_stimOn)), 'Color', ones(1,3).*0.5)

set(fig_grand, 'Color', 'w', 'Position', [100   100   1594   924])
setCell_str = sprintf('%d_', setCell);
setCell_str = setCell_str(1:end-1);
print(fig_grand, fullfile(dirFig, sprintf('%s_FOV%d_session%d_exCells_%s_MOV_IMG_RS', ...
    nameSubj, FOV_ID, iSession, setCell_str)), '-r0', '-depsc');



%%% 2. draw IMG, MOV, RS responses separately: 3 or 4 neurons for each  BPM>DFL, BPM<DFL, BPM=DFL 
% set of the example cells 
setCell = [40 9 73 21 4 35 79 36 23]; %[40 9 73 82 21 4 35 41 79 36 23 29]; %[40 9 73 82 21 4 35 41 79 36 23 29];

% prepare the cell traces
matTS_BPM = [];
for iCell = 1:length(setCell)
    matTS_BPM(:, iCell) = cat(1, tS_run(1).tS_trial(setCell(iCell),:).matTS); % this is C_raw
end

iRun = 1; iMovie = 1;
matTS_DFL = [];
matTS_DFL = tS_session(iMovie).matTS_C_raw(:, setCell, iRun);

load(fullfile(dirProcdata_session, 'RS_ts.mat'))
matTS_RS = [];
matTS_RS = tSeries_RS(1).C_raw(setCell, :)';


cMap_lineColor = [102 194 165; 252 141 98; 141 106 203]./255;
cMap_lineColor2 = [27 158 119; 217 95 2; 117 112 179]./255;

% plot BPM responses of a set of example cells with stim ON marks
fig_BPM = figure;
p = plot(matTS_BPM+repmat([1:length(setCell)].*3, size(matTS_BPM, 1), 1), 'LineWidth', 2);
axis tight
line(repmat(11:46:2208, 2, 1), repmat(get(gca, 'ylim')', 1, 48), 'Color', ones(1,3).*0.5)
xlim([1 1200])
set(p(1:3), 'Color', cMap_lineColor(1,:))
set(p(4:6), 'Color', cMap_lineColor(2,:))
set(p(7:9), 'Color', cMap_lineColor(3,:))
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'k', 'LineWidth', 2)
set(fig_BPM, 'Color', 'w', 'Position', [385   185   500   940])
% print(fig_BPM, fullfile(dirFig, sprintf('%s_FOV%d_session%d_exCells_BPMresp_group', ...
%     nameSubj, FOV_ID, iSession)), '-r0', '-depsc');


% plot DFL responses of a set of example cells
fig_DFL = figure;
p = plot(matTS_DFL+repmat([1:length(setCell)].*3, size(matTS_DFL, 1), 1), 'LineWidth', 2);
axis tight
line(repmat(200:200:1200, 2, 1), repmat(get(gca, 'ylim')', 1, 6), 'Color', ones(1,3).*0.5)
xlim([1 1200])
set(p(1:3), 'Color', cMap_lineColor(1,:))
set(p(4:6), 'Color', cMap_lineColor(2,:))
set(p(7:9), 'Color', cMap_lineColor(3,:))
hold on;
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'k', 'LineWidth', 2)
set(fig_DFL, 'Color', 'w', 'Position', [100   100   500   940]) 
% print(fig_DFL, fullfile(dirFig, sprintf('%s_FOV%d_session%d_exCells_DFLresp_group', ...
%     nameSubj, FOV_ID, iSession)), '-r0', '-depsc');


% plot RS responses of a set of example cells
fig_RS = figure;
p = plot(matTS_RS+repmat([1:length(setCell)].*3, size(matTS_RS, 1), 1), 'LineWidth', 2);
axis tight
xlim([2100 2100+1199])
set(p(1:3), 'Color', cMap_lineColor(1,:))
set(p(4:6), 'Color', cMap_lineColor(2,:))
set(p(7:9), 'Color', cMap_lineColor(3,:))
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'LineWidth', 2)
set(fig_RS, 'Color', 'w', 'Position', [100   100   500   940])
% print(fig_RS, fullfile(dirFig, sprintf('%s_FOV%d_session%d_exCells_RSresp_group', ...
%     nameSubj, FOV_ID, iSession)), '-r0', '-depsc');


%% Mark these example cells in the FOV
figure;
subplot('Position', [0 0 1 1]);
image(ones(size(neuron.Cn)).*255);
colormap(gray);
hold on

thr = 0.5; %0.2;
cellColor = ones(1,3).*0.75; %'k'; %cMap_line(iSession, :); %'m'; %'c'; 
widthContour = 2;
[d1,d2] = size(neuron.Cn);
indCellValid_session = cellIDAcrossDay(~isnan(cellIDAcrossDay(:,iSession)), iSession);

CC = cell(length(indCellValid_session),1);
CR = cell(length(indCellValid_session),2);
for iCC = 1:length(indCellValid_session)
    i = indCellValid_session(iCC); %sortedIndCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
truesize;

hold on;
cmap_cell = repmat(colororder, 2, 1); %colororder; % cmap_cell = cool(size(neuron.A, 2));
[a, setIDCell] = ismember(setCell, indCellValid_session); %ismember(ind(1:5), indCellValid_session);

for iCC = 1: length(setIDCell)
    i = setIDCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', cmap_cell(iCC,:), 'linewidth', widthContour+0.5);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end

line([1 33], [265 265], 'Color', 'k', 'LineWidth', 5)

setCell_str = sprintf('%d_', setCell);
setCell_str = setCell_str(1:end-1);
print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_gray_scaleBar_exCells_%s', ...
    nameSubj, FOV_ID, thr*10, iSession, setCell_str)), '-r0', '-depsc');



%% Previous: Example time courses
% dateSession = setDateSession{1};
% dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% dirPreproc = fullfile(dirProcdata_session, '_preproc');
% 
% load(fullfile(dirProcdata_session, 'DFL_ts.mat'))
% 
% % select the cells
% pnrs_dfl = max(tSeries_DFL(1).C, [], 2)./std(tSeries_DFL(1).C_raw-tSeries_DFL(1).C, 0, 2); % peak amplitude divided by noise std
% [a, ind] = sort(pnrs_dfl, 'descend');
% % snrs = var(tSeries_DFL(1).C, 0, 2)./var(tSeries_DFL(1).C_raw-tSeries_DFL(1).C, 0, 2);
% % [a, ind] = sort(snrs, 'descend');
% 
% fig_movie = figure;
% set(fig_movie, 'Color', 'w')
% tlen = 1200;
% plot(tSeries_DFL(1).C_raw(ind(1:5), 1:tlen)'+repmat([1:5].*5, tlen, 1), 'LineWidth', 2)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
% 
% print(fig_movie, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_DFL1_DFLsnr5', dateSession)), '-depsc')
% 
% 
% load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
% 
% pnrs_bpm = max(tSeries_BPM(1).C, [], 2)./std(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2); % peak amplitude divided by noise std
% [a, ind] = sort(pnrs_bpm, 'descend');
% % snrs = var(tSeries_BPM(1).C, 0, 2)./var(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2);
% % [a, ind] = sort(snrs, 'descend');
% 
% fig_BPM = figure;
% set(fig_BPM, 'Color', 'w')
% tlen = 1200;
% plot(tSeries_BPM(1).C_raw(ind(1:5), 21:20+tlen)'+repmat([1:5].*5, tlen, 1), 'LineWidth', 2)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
% 
% print(fig_BPM, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_BPM1_DFLsnr5', dateSession)), '-depsc')







