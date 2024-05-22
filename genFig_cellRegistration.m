% genFig_cellRegistration.m
%
% 2024/05/20 SHP
% visualize the longitudinal registration results
%   - a histogram of tracking duration?


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

iSubj = 2; %1; %1; %2; %1; %2; %1;

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

% % aligned cells movie TS and spatial info
% fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
% load(fname_caTSFOV, 'cellTS', 'cellPix')


%% Plot the histogram of cells detected across sessions with cumulative distribution
flagCell = ~isnan(cellIDAcrossDay);

figure;
yyaxis left
h = histogram(sum(flagCell, 2));
yyaxis right
plot(cumsum(h.Values), 'o-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w')
set(h, 'lineWidth', 2)
set(gca, 'TickDir', 'out')
set(gca, 'XTick', 1:12)
box off
set(gcf, 'Color', 'w')
set(gcf, 'Position', [100 100 500 500])
set(gca, 'FontSize', 15, 'LineWidth', 2)
% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_cellRegistration_histogram', nameSubj, FOV_ID)), '-r300', '-depsc')


%% Plot the cell registration results over days, not number of sessions
setDateSession_datenum = datenum(setDateSession, 'yyyymmdd');
nDaysBetweenSessions = diff(setDateSession_datenum);

matDaysAcRegistration = [];
nSessionRegistered = sum(flagCell, 2);
for iCell = 1:length(nSessionRegistered)
    if nSessionRegistered(iCell) > 1
        locSession = find(flagCell(iCell, :)>0);
        nDays = setDateSession_datenum(locSession(end)) - setDateSession_datenum(locSession(1));
        matDaysAcRegistration(iCell,1) = nDays;
    else
        matDaysAcRegistration(iCell,1) = 0;
    end
end


figure;
set(gcf, 'Color', 'w', 'Position', [100 100 500 500])
hh = histogram(matDaysAcRegistration(matDaysAcRegistration>0));
set(gca, 'TickDir', 'out', 'Box', 'off')
set(gca, 'FontSize', 15, 'LineWidth', 2)




