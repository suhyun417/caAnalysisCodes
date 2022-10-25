% genFig_2022SFN.m
%
% 2022/09/26 SHP
%   - generate figures for 2022 SFN talk
%   - part of this script is to generate a figure for 2022 NIMH Training Day
%  poster
% 

clear all;

%% Settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
%     addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/NIFVAULT/projects';
        directory.procdata = '/Volumes/NIFVAULT/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/NIFVAULT/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/projects';
        directory.procdata = '/nifvault/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/nifvault/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end
dirFig = '/nifvault/projects/parksh/0Marmoset/Ca/_labNote/_figs';

%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

nameSubj = 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = 1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%%


%% Read source data
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


% get the contours and image field of view
% neuron_b = neuron.batches{1}.neuron;

% Generate cell location map within FOV
thr = 0.5; % the lower the smaller (more centralized) the contour
cellColor = [1 1 1];
widthContour = 1;
[d1,d2] = size(neuron.Cn);

figure;
imagesc(zeros(d1, d2)); % background
colormap(gray);
caxis([0 0.1]);
hold on;

CC = cell(size(neuron.A, 2),1);
CR = cell(size(neuron.A, 2),2);
for i = 1:size(neuron.A ,2)
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
title(sprintf('%s: %s', nameSubj, dateSession))
% %save
% print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');

% end
%
Coor = neuron.get_contours(thr); 
CC = Coor;
    for i = 1:size(Aor,2)
%         cont = medfilt1(Coor{i}')';
        cont = Coor{i}; 
        if size(cont,2) > 1
            plot(cont(1,1:end),cont(2,1:end),'Color',cmap(i+size(Aor,2),:), 'linewidth', ln_wd); hold on;
        end
    end

% %% 2022 Training Day figures
% %% Intro figure: three different movie-driven ts of macaque face cells
% load('/nifvault/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'matTS_FP')
% 
% setCellIDs = {'06Dav', '25Dav', '27Dav'}; % three example neurons from the face patch ML
% tLoc = find(contains(matTS_FP.catChanID, setCellIDs));  
% 
% cellTS = matTS_FP.matFR_SU_10hz(1:3000, tLoc);
% cellTS_rs = resample(cellTS, 1, 10);
% 
% % Plot
% cMap_cell = [96 216 54; 255 67 161; 0 162 253]./255; % three example neurons
% 
% fig_intro = figure;
% set(fig_intro, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1000 500 650 480]);
% % taxis = 2.4:2.4:900;
% for iCell = 1:3
%     figure(fig_intro);
%     plot(cellTS_rs(:, iCell)+1*(iCell-1), 'LineWidth', 2, 'Color', cMap_cell(iCell,:));
%     hold on;
% end
% set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
% print(fig_intro, fullfile(dirFig, sprintf('macaqueCellMovie1_%s_%s_%s', setCellIDs{1}, setCellIDs{2}, setCellIDs{3})), ...
%     '-depsc')
% 
% 
% %% Correlation values
% matR = corr(cellTS, fmriTS, 'rows', 'complete', 'type', 'spearman');
% % matR = matR*(-1); %becuase of the MION
% 
% %% face regressor ts
% taxis = 2.4:2.4:900;
% 
% faceRegressor_mion([1:7,126:126+7, 251:251+7]) = NaN;
% faceRegressor_mion_norm = (faceRegressor_mion-nanmean(faceRegressor_mion))./nanstd(faceRegressor_mion);
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 100 650 160]);
% plot(taxis, faceRegressor_mion_norm, '-', 'LineWidth', 1, 'Color', [154 205 50]./255);
% set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
% set(gca, 'YColor', 'none')
% print(gcf, fullfile(dirFig, 'multipleFP_sFig_faceRegressor_mion'), '-depsc')

