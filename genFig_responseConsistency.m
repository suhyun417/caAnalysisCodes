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
flagCell = ~isnan(cellIDAcrossDay);

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


%%
locWeek = find(matDaysAcRegistration>=14);
length(locWeek)

% gather averaged TS of 1st session and last session that is at least 2
% weeks apart
matAvgTS_1 = []; matAvgTS_2 = []; tempInfoCell = [];
for iCell = 1:length(locWeek)
    idCell = locWeek(iCell);
    locSession = find(flagCell(idCell, :)>0);
    nDays = setDateSession_datenum(locSession(end)) - setDateSession_datenum(locSession(1));

    setT = []; 
    setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);

    matAvgTS_1(:,iCell) = mean(cellTS(idCell).matTS_movie1(setT(1,1):setT(1,2),:),1)';
    matAvgTS_2(:,iCell) = mean(cellTS(idCell).matTS_movie1(setT(end,1):setT(end,2),:),1)';
    tempInfoCell(iCell, 1) = idCell;
    tempInfoCell(iCell, 2) = nDays;

end

[matR] = corr(matAvgTS_1, matAvgTS_2, 'type', 'Spearman');
setR = diag(matR);
figure
hR = histogram(setR);
figure
plot(tempInfoCell(:,2), setR, 'o')

% quick test
for iCell = 1:length(locWeek)
    figure(3); cla;
    plot(matAvgTS_1(:,iCell), 'ro-');
    hold on
    plot(matAvgTS_2(:,iCell), 'bo-');
    input('')
end


%%
figtemp = figure;
set(gcf, 'color', 'w', 'Position', [401    35   378   862])

tempSetC = find(matDaysAcRegistration>25);
for iC = 1:length(tempSetC)

idCell = tempSetC(iC); %34; %23;

setT = []; avgT=[];
setT = cat(2, cat(1, 1, cellTS(idCell).nTrial1_set(1:end-1)+1), cellTS(idCell).nTrial1_set);
for iSS = 1:length(setT)
    avgT(:,iSS) = mean(cellTS(idCell).matTS_movie1(setT(iSS,1):setT(iSS,2),:),1)';
end


stringNameSession= cat(1, setDateSession{cellTS(idCell).idAcrossSession(:,1)}); % get date
cMap = cool(length(setT));

% figtemp = figure;
% set(gcf, 'color', 'w', 'Position', [401    35   378   862])
figure(figtemp); cla;
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
title(sprintf('%s: Cell ID %d', nameSubj, idCell))
% print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_CellID%d_avgTS_eachSession',nameSubj, FOV_ID, idCell)), '-depsc')

input('')
end


%% In case I need the histogram values
% figure;
% h = histogram(sum(flagCell, 2));
% 
% plot(cumsum(fliplr(h.Values)), 'o-')
% hold on
% text((1:length(h.Values))'-0.25, (cumsum(fliplr(h.Values))+10)', num2str(cumsum(fliplr(h.Values))'))
% set(gca, 'XTick', 1:length(h.Values), 'XTickLabel', length(h.Values):-1:1)
% xlabel('Number of detected sessions')
% ylabel('Cumulative number of cells')
% set(gca, 'TickDir', 'out', 'Box', 'off')
% set(gca, 'FontSize', 15, 'LineWidth', 2)
% % print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_cellRegistration_cumNumOverSessions', nameSubj, FOV_ID)), '-r300', '-depsc')



