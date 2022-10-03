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

%% 2022 Training Day figures
%% Intro figure: three different movie-driven ts of macaque face cells
load('/nifvault/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'matTS_FP')

setCellIDs = {'06Dav', '25Dav', '27Dav'}; % three example neurons from the face patch ML
tLoc = find(contains(matTS_FP.catChanID, setCellIDs));  

cellTS = matTS_FP.matFR_SU_10hz(1:3000, tLoc);
cellTS_rs = resample(cellTS, 1, 10);

% Plot
cMap_cell = [96 216 54; 255 67 161; 0 162 253]./255; % three example neurons

fig_intro = figure;
set(fig_intro, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1000 500 650 480]);
% taxis = 2.4:2.4:900;
for iCell = 1:3
    figure(fig_intro);
    plot(cellTS_rs(:, iCell)+1*(iCell-1), 'LineWidth', 2, 'Color', cMap_cell(iCell,:));
    hold on;
end
set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
print(fig_intro, fullfile(dirFig, sprintf('macaqueCellMovie1_%s_%s_%s', setCellIDs{1}, setCellIDs{2}, setCellIDs{3})), ...
    '-depsc')


%% Correlation values
matR = corr(cellTS, fmriTS, 'rows', 'complete', 'type', 'spearman');
% matR = matR*(-1); %becuase of the MION

%% face regressor ts
taxis = 2.4:2.4:900;

faceRegressor_mion([1:7,126:126+7, 251:251+7]) = NaN;
faceRegressor_mion_norm = (faceRegressor_mion-nanmean(faceRegressor_mion))./nanstd(faceRegressor_mion);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 100 650 160]);
plot(taxis, faceRegressor_mion_norm, '-', 'LineWidth', 1, 'Color', [154 205 50]./255);
set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
set(gca, 'YColor', 'none')
print(gcf, fullfile(dirFig, 'multipleFP_sFig_faceRegressor_mion'), '-depsc')

