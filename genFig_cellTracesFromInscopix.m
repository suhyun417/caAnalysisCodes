% genFig_cellTracesFromInscopix.m
%

dirDataHome = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/20191010'; 
% dirDataHome = '/Users/parks20/Documents/Inscopix_Projects/20191010_Max_timeseries/';
% fname = '20191010_Max_124810_CellTraces.csv'; 
fname = '20191010_Max_125432_CellTraces.csv'; %'20191010_Max_124810_CellTraces.csv';
fname_props = '20191010_Max_125432_CellTraces-props.csv'; %'20191010_Max_124810_CellTraces.csv';


% dirDataHome = '/Users/parks20/Documents/Inscopix_Projects/20191004_Max/';
% fname = '20191004_Max_133333_CellTraces.csv'; 

T = readtable(fullfile(dirDataHome, fname));
status = T{1,:};
tf = contains(status, 'accepted');
locValidCells = find(tf>0);
cellNames = T.Properties.VariableNames(locValidCells);

tempVals = csvread(fullfile(dirDataHome, fname), 2, 0);
tAxis = tempVals(:,1); % time axis
matTS = tempVals(:, locValidCells); %T{2:end, locValidCells};

nCells = size(matTS, 2);

% iCell = 1;
% 
% figure;
% plot(tAxis, matTS(:, iCell), '-', 'LineWidth', 2);
% title(sprintf('Filename: %s, Cell Name: %s', fname, cellNames{iCell}))


% dirBHV = '/Users/parks20/ResearchProjects/0MARMOSET/0FaceNeurons/_data/invivoCalciumImaging/Max/20191004_Max/';
% setFileName = {'191004_133325_basicPreferenceMapping_userloop_Max.bhv2', ...
% '191004_135020_basicPreferenceMapping_userloop_Max.bhv2', ...
% '191004_135936_basicPreferenceMapping_userloop_Max.bhv2', ...
% '191004_141730_basicPreferenceMapping_userloop_Max.bhv2', ...
% '191004_142434_basicPreferenceMapping_userloop_Max.bhv2'};

dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %'/Users/parks20/ResearchProjects/0MARMOSET/0FaceNeurons/_data/invivoCalciumImaging/Max/20191010_Max'; % 20191004_Max/';
setFileName = {'191010_124747_basicPreferenceMapping_userloop_Max.bhv2', ...
'191010_125402_basicPreferenceMapping_userloop_Max.bhv2',...
'191010_130246_basicPreferenceMapping_userloop_Max.bhv2',...
'191010_130811_basicPreferenceMapping_userloop_Max.bhv2',...
'191010_132228_basicPreferenceMapping_userloop_Max.bhv2',...
'191010_132323_basicPreferenceMapping_userloop_Max.bhv2'};

% for iF = 1:length(setFileName)
iF = 2;
    [data,MLConfig,TrialRecord] = mlconcatenate_woUserVars(fullfile(dirBHV, setFileName{iF}));
    t_sendTTL = data.BehavioralCodes(1).CodeTimes(2);
    t_startCa = find(data.AnalogData.Button.Btn1>0, 1);
    dT = t_startCa-t_sendTTL;
    fprintf(1, '         Filename: %s, t_sendTTL: %3.3f, t_startCa: %3.3f, dT: %3.3f\n', setFileName{iF}, t_sendTTL, t_startCa, dT);
% end


t_endTTL = data.BehavioralCodes(40).CodeTimes(data.BehavioralCodes(40).CodeNumbers==990);
t_trialStart_abs = data.AbsoluteTrialStartTime;
catCodeTimes = cat(1, data.BehavioralCodes.CodeTimes);
catCodeNumbers = cat(1, data.BehavioralCodes.CodeNumbers);
t_stimOnset = catCodeTimes(catCodeNumbers==20);
t_blankOnset = catCodeTimes(catCodeNumbers==15);
t_afterOff = catCodeTimes(catCodeNumbers==25)+3000;

t_onset = t_stimOnset; %t_blankOnset;
t_end = t_afterOff;

fs = 15; %10;
delay = 0; %0.2; %t_startCa/1000; %0.2;
t_onset_adj = t_onset/1000-delay;
locStart = floor(t_onset_adj./(1/fs));
locEnd = floor((t_end/1000-delay)./(1/fs));  %cat(1, locOnset(2:end)-1, floor((t_endTTL/1000-delay)./(1/fs)));

clear matTS_trial
for iTrial = 1:length(locStart)
    matTS_trial{iTrial} = matTS(locStart(iTrial):locStart(iTrial)+49, :); %locEnd(iTrial), :); %locOnset_end(iTrial), :);
end

cond = cat(1, data.Condition(data.TrialError>1));
[idCond, locTrial] = sort(cond);

setCond = unique(idCond);
nCond = length(setCond);

% clear matTS_sorted %= [];
matTS_sorted = [];
for iCond = 1:nCond
    curLoc = locTrial(idCond==setCond(iCond));
    catCond = cat(3, matTS_trial{curLoc});
    matTS_sorted = cat(4, matTS_sorted, catCond); % time x cell x trials x conditions
end

setCellInd = [16 10 19];
cmap = [255 47 146; 142 250 0; 118 214 255]./255;
mergeCond(1).cond = [1 2]; % faces: hf & mf
mergeCond(1).condname = 'human face, marmoset face';
mergeCond(2).cond = [4 6]; % objs: scene & objects
mergeCond(2).condname = 'scene, familiar objects';
mergeCond(3).cond= [9 10]; % motion: gratings & rdm
mergeCond(3).condname = 'grating, random dot motion';


for iCell = 1: length(setCellInd)
    
    for iMC = 1:length(mergeCond)
        tempAll=[];
        tempAll = matTS_sorted(:, setCellInd(iCell), :, ismember(setCond, mergeCond(iMC).cond));
        
        exampleCells(iCell, iMC).matTS = reshape(squeeze(tempAll), size(tempAll, 1), size(tempAll, 3)*size(tempAll, 4));
        exampleCells(iCell, iMC).meanTS = mean(exampleCells(iCell, iMC).matTS, 2);
        exampleCells(iCell, iMC).indCell = setCellInd(iCell);
        exampleCells(iCell, iMC).cellID = cellNames{setCellInd(iCell)};
        exampleCells(iCell, iMC).mergeCond = mergeCond(iMC).cond;
        exampleCells(iCell, iMC).mergeCond_name = mergeCond(iMC).condname;
    end
end

% save figures
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';
for iCell = 1:length(setCellInd)
    tt = strsplit(fname, '_');
    for iMC = 1:length(mergeCond)
        figure_ex_trials = figure;
        set(figure_ex_trials, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 140 300])
        plot(exampleCells(iCell, iMC).matTS, 'LineWidth', 2, 'Color', cmap(iCell,:))
        ylim([-6 12])
        set(gca, 'Box', 'off', 'TickDir', 'out')
%         title(sprintf('ExampleCaTraces_%s_%s_mergeCond%d_matTS', strjoin(tt(1:3), '_'), exampleCells(iCell, iMC).cellID, iMC))
        
        figure_ex_mean = figure;
        set(figure_ex_mean, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 140 300])
        plot(exampleCells(iCell, iMC).meanTS, 'LineWidth', 2, 'Color', cmap(iCell,:))
        ylim([-1 2])
        set(gca, 'Box', 'off', 'TickDir', 'out')
%         title(sprintf('ExampleCaTraces_%s_%s_mergeCond%d_meanTS', strjoin(tt(1:3), '_'), exampleCells(iCell, iMC).cellID, iMC))
        
        % save
        print(figure_ex_trials, fullfile(dirFig, sprintf('ExampleCaTraces_%s_%s_mergeCond%d_matTS', strjoin(tt(1:3), '_'), exampleCells(iCell, iMC).cellID, iMC)), '-depsc');
        print(figure_ex_mean, fullfile(dirFig, sprintf('ExampleCaTraces_%s_%s_mergeCond%d_meanTS', strjoin(tt(1:3), '_'), exampleCells(iCell, iMC).cellID, iMC)), '-depsc');
    end
end

for iCell = 1:nCells
    % curTS_trial = matTS_sorted(:, iCell, :, iCond);
    
    cmap = jet(nCond);
    fig_condition=figure;
    for iCond = 1:nCond
        curTS_trial = matTS_sorted(:, iCell, :, iCond);
        figure(fig_condition);
        subplot(nCond, 1, iCond);
        plot(squeeze(curTS_trial), 'Color', cmap(iCond, :), 'LineWidth', 2);
        ylim([-1 20])
        if iCond == 1
            title(sprintf('Filename: %s, Cell Name: %s', fname, cellNames{iCell}))
        end
    end
    fprintf(1, 'Filename: %s, Cell Name: %s\n', fname, cellNames{iCell})
    input('')
end


% movie watching
dirDataHome = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/20191010'; 
% dirDataHome = '/Users/parks20/Documents/Inscopix_Projects/20191010_Max_timeseries/';
% fname = '20191010_Max_124810_CellTraces.csv'; 
fname_movie = '20191010_Max_131421_CellTraces.csv'; %'20191010_Max_124810_CellTraces.csv';

% dirDataHome = '/Users/parks20/Documents/Inscopix_Projects/20191004_Max/';
% fname = '20191004_Max_133333_CellTraces.csv'; 

T_movie = readtable(fullfile(dirDataHome, fname_movie));
status = T_movie{1,:};
tf = contains(status, 'accepted');
locValidCells = find(tf>0);
cellNames = T_movie.Properties.VariableNames(locValidCells);

tempVals = csvread(fullfile(dirDataHome, fname_movie), 2, 0);
tAxis = tempVals(:,1); % time axis
matTS_movie = tempVals(:, locValidCells); %T{2:end, locValidCells};

nCells = size(matTS_movie, 2);

setCellInd_movie = [14 20 7];
figure_movie = figure;
set(figure_movie, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 900 350])
for iCell = 1:length(setCellInd_movie)
    SP(iCell) = subplot(length(setCellInd_movie), 1, iCell);
    plot(matTS_movie(:, setCellInd_movie(iCell)), 'LineWidth', 2, 'Color', cmap(iCell,:));
end
% ylim([-6 12])
set(SP, 'Box', 'off', 'TickDir', 'out', 'YLim', [-6 12], 'XLim', [0 1800])

tt = strsplit(fname_movie, '_');
print(figure_movie, fullfile(dirFig, sprintf('ExampleCaTraces_movie_%s_%d%d%d', strjoin(tt(1:3), '_'), setCellInd_movie)), '-depsc');

