% genFig_responseConsistency.m
%
% 2024/05/13 SHP
% visualize the consistency across trials/sessions
% focus on finding the most effective way to show consistency of
%   - example neurons, populations
%   - over different timescale
%   - over different aspects, e.g. selectivity, fluctuations etc.


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

% aligned cells movie TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')


%%
idCell = 34; %23;

setT = []; avgT=[];
setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);
for iSS = 1:length(setT)
    avgT(:,iSS) = mean(cellTS(idCell).matTS_movie1(setT(iSS,1):setT(iSS,2),:),1)';
end


stringNameSession= cat(1, setDateSession{cellTS(idCell).idAcrossSession(:,1)}); % get date
cMap = cool(length(setT));

figtemp = figure;
set(gcf, 'color', 'w', 'Position', [401    35   378   862])
for iSS = 1:length(setT)
figure(figtemp);
ydata = avgT(:,iSS)+50*(iSS-1);
plot(ydata, 'color', cMap(iSS,:));
ytick(iSS) = mean(ydata);
hold on;
end
axis tight
set(gca, 'TickDir', 'out', 'box', 'off')
set(gca, 'YTick', ytick, 'YTicklabel', stringNameSession)
set(gca,'XTick', 200:200:1200, 'XTickLabel', 20:20:120)
% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession',nameSubj, FOV_ID, idCell)), '-depsc')




