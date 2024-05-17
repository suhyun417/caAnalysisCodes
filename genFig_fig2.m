% genFig_fig2.m
% 
% 2024/04/22 SHP
% Generate Figure 2 of marmoset calcium imaging manuscript
% a) example FOV and calcium traces
% b) longitudinal registration of the FOVs in two animals: extracted
% sources from each day, superimposed source boundaries over days
% c) entire population of cells (colored to indicate longitudinal
% tracking?)
% d) Histogram of tracking durations.  Need to think about whether total days, one for consecutive, or longest span.


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

% aligned cells TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')

%% A. Example FOV and representative time series
%%% FOV
iSession = 1; %:length(setDateSession) %;
dateSession = setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

% dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

load(sprintf('%s/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', directory.dirProcdata, nameSubj, dateSession))
% load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data
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

print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_gray_scaleBar_r0', ...
    nameSubj, FOV_ID, thr*10, iSession)), '-r0', '-depsc');


%%% example time courses
dateSession = setDateSession{1};
dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

load(fullfile(dirProcdata_session, 'DFL_ts.mat'))

% pnrs = max(tSeries_DFL(1).C, [], 2)./std(tSeries_DFL(1).C_raw-tSeries_DFL(1).C, 0, 2); % peak amplitude divided by noise std
% [a, ind] = sort(pnrs, 'descend');
snrs = var(tSeries_DFL(1).C, 0, 2)./var(tSeries_DFL(1).C_raw-tSeries_DFL(1).C, 0, 2);
[a, ind] = sort(snrs, 'descend');

fig_movie = figure;
set(fig_movie, 'Color', 'w')
tlen = 1200;
plot(tSeries_DFL(1).C_raw(ind(1:5), 1:tlen)'+repmat([1:5].*5, tlen, 1), 'LineWidth', 2)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
print(fig_movie, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_DFL1_DFLsnr5', dateSession)), '-depsc')


load(fullfile(dirProcdata_session, 'BPM_ts.mat'))

% pnrs = max(tSeries_BPM(1).C, [], 2)./std(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2); % peak amplitude divided by noise std
% [a, ind] = sort(pnrs, 'descend');
% snrs = var(tSeries_BPM(1).C, 0, 2)./var(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2);
% [a, ind] = sort(snrs, 'descend');

fig_BPM = figure;
set(fig_BPM, 'Color', 'w')
tlen = 1200;
plot(tSeries_BPM(1).C_raw(ind(1:5), 21:20+tlen)'+repmat([1:5].*5, tlen, 1), 'LineWidth', 2)
set(gca, 'Box', 'off', 'TickDir', 'out', 'YColor', 'w', 'XTick', 0:200:tlen, 'XTickLabel', 0:20:tlen/10, 'LineWidth', 2)
print(fig_BPM, fullfile(dirFig, sprintf('exampleTS_Tabla_%s_BPM1_DFLsnr5', dateSession)), '-depsc')


%%% Mark these example cells in the FOV
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
cmap_cell = colororder; % cmap_cell = cool(size(neuron.A, 2));
[a, setIDCell] = ismember(ind(1:5), indCellValid_session);

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
print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_gray_scaleBar_exCells', ...
    nameSubj, FOV_ID, thr*10, iSession)), '-r0', '-depsc');


%% B & C. longitudinal registration (cells from all the sessions) & its aligned version
% nameSubj = 'Tabla';
% dateSession = '20191125'; % '20191113';

cMap_line = cool(length(setDateSession));
% fig_contours = figure;
for iSession = 1:length(setDateSession) %;
dateSession = setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

% dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

load(sprintf('%s/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', directory.dirProcdata, nameSubj, dateSession))
% load(sprintf('/nifvault/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data
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

print(fullfile(dirFig, sprintf('%s_FOV%d_validCellContour_thr0p%d_session%d_gray_scaleBar', ...
    nameSubj, FOV_ID, thr*10, iSession)), '-depsc');

end



%% D. FOV of entire population of cells
%%% FOV of all the cells

% tempA = cat(2, cellPix(:).repPix);
% tempA(~isnan(tempA)) = 1;
% 
% % % selected cells
% % tempA = cat(2, cellPix(indCellValid).repPix);
% % tempA(~isnan(tempA)) = 1;
% % % tempA(:, indCellValid) = tempA(:, indCellValid).*10;
% 
% imgCells = sum(tempA, 2, 'omitnan');
% imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));
% 
% figure;
% set(gcf, 'Color', 'w')
% subplot('Position', [0 0 1 1])
% imagesc(imgCells_2d)
% colormap(gray)
% set(gca, 'CLim', [0 0.5])
% truesize
% box off
% axis off
% % print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_cellMap_allCells_truesize', nameSubj, FOV_ID)), '-depsc')


%% D. 











