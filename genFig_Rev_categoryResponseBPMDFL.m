% genFig_Rev_categoryResponseBPMDFL.m
%
% For Imaging Neuroscience revision
% 2025/11/19 modified by SHP
%   - this code now deals with the plotting categorical responses from BPM runs 
%   for the revised manuscript
%   - the part doing normalization & save the matrix is now separated to
%   "runSaveCaTS_BPMsorted_baselineNormalization.m"
% 2025/11/13 created by SHP
%   - create bar graph showing categorical responses to images for example neurons
%   - possibly scatter plot to compare image vs. movies
%   - For that, stimulus response needs to be recomputed using baseline subtraction
%   - so instead of using currently saved normalized amplitude and
%   baseline, it needs to be re-computed

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

iSubj = 1; %2; %1; %2; %1; %2; %1; %2; %1;

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
fname_caTSFOV_DFL = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_DFL, 'cellTS', 'cellPix')
cellTS_DFL = cellTS;
clear cellTS

fname_caTSFOV_BPM = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_BPMsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV_BPM, 'cellTS')
cellTS_BPM = cellTS;
clear cellTS

% fname_caTSFOV_RS = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
% load(fname_caTSFOV_RS, 'cellTS_RS') %, 'resultsRS')


%% gathering the normalized stimulus responses to images and categories
indCellValid = find(cat(1, cellTS_DFL.nTrial1_total)>8); %

clear matAvg*
for iC = 1:length(indCellValid)
    
    clear resp* matTS_norm* tempCat
    iCell = indCellValid(iC);
    
    tempCat = cat(1, cellTS_BPM(iCell,:).matTS_baselinenorm_avg);
    matTS_bnorm_cat = reshape(tempCat', [46 5 5]);
    matTS_bnorm_cat_avg = squeeze(mean(matTS_bnorm_cat, 2));
    
    resp_avg_img = mean(tempCat(:, 20:30), 2)';
    resp_avg_cat = mean(matTS_bnorm_cat_avg(20:30, :));
    

    matAvg_img(iC, :) = resp_avg_img; % cell x image
    matAvg_cat(iC, :) = resp_avg_cat; % cell x category
end

figure
imagesc(matAvg_cat)
colormap(hot)
colorbar
set(gca, 'CLim', [-1 1]*1.96)
figure
imagesc(matAvg_img)
colormap(hot)
colorbar
set(gca, 'CLim', [-1 1]*1.96)

% selected cells for Figure 2 example traces: original cell IDs are from
% Tabla Session 1
setCell_org = [9 73 82 29 36 79 41 4 35 21]; 
[a, setCell_ind] = ismember(setCell_org, cellIDAcrossDay(:,1));

matAvg_cat_setCell = matAvg_cat(setCell_ind, :);
% L = legend('HF', 'MF', 'NFO', 'FO', 'RDM');
figure;
tiledlayout(10,1)
for iB = 1:10
    ax(iB) = nexttile;
    bar(ax(iB), matAvg_cat_setCell(iB,:));
end

% For selected cells, let's recollect the image responses
for iSC = 1:length(setCell_ind)

    iCell = setCell_ind(iSC);

    tempCat = cat(1, cellTS_BPM(iCell,:).matTS_norm_avg);
    matTS_bnorm_cat = reshape(tempCat', [46 5 5]);
    matTS_bnorm_cat_avg = squeeze(mean(matTS_bnorm_cat, 2));
    
    resp_avg_img = mean(tempCat(:, 20:30), 2)';
    resp_avg_cat = mean(matTS_bnorm_cat_avg(20:30, :));





















