% analCa_ttestBPM.m
%
% script to evaluate visually responsive neurons based on responses 
% from flashing image (BPM) runs using t-test with Bonferroni correction
% 2020/03/11 SHP

clear all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

    
nameSubj ='Tabla'; % setNameSubj{iSubj}; %'Max'; %

% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[cc, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = cc(2:end); % 1st one is always empty
nSession = length(setDateSession);

close all;

iSession = 1; %12; %1;
dateSession = setDateSession{iSession};
% dateSession = '20191113'; %setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% %% Read source data and compute center coordinates of cells
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
% [center] = neuron.estCenter();
% center = center(validIndCell, :);

%% t-test to compare visual responses and baseline
load(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession)), 'tS_session_stim')

a = 0.05;
for iCell = 1:size(tS_session_stim, 1)
    for iStim = 1:size(tS_session_stim, 2)
        [h, p] = ttest(tS_session_stim(iCell, iStim).matAvgAmp_norm, tS_session_stim(iCell, iStim).matAvgAmp_b_norm, 'alpha', a/size(tS_session_stim,2));
        tval(iStim, iCell) = h;
        pval(iStim, iCell) = p;
    end
end

find(sum(tval)<1)

% Columns 1 through 15
% 
%     11    12    13    14    16    17    18    20    22    29    32    38    44    52    53
% 
%   Columns 16 through 26
% 
%     58    70    72    77    83    85    98    99   107   109   116

%   Columns 1 through 15
% 
%      7    13    14    17    18    32    33    38    46    58    71    72    83    85    99
% 
%   Columns 16 through 17
% 
%    107   109
