% analCa_coherogram.m

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

iSession = 1;
dateSession = setDateSession{iSession};
% dateSession = '20191113'; %setDateSession{iSession};

dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');

% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession))
% load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession))

%% Read source data and compute center coordinates of cells
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

[center] = neuron.estCenter();
center = center(validIndCell, :);

%% Pairwise coherogram
load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_run')

params.pad = 1;
params.Fs = 10;
params.trialave = 1;
% params.err      = [0 0.05];

idCell1 = 105;
idCell2 = 108;
data1 = []; data2=[];
for iRun = 1:length(tS_run)
curS = tS_run(iRun).tS_trial;
curTS1 = cat(1, curS(idCell1, 1:minNumTrial).matTS);
curTS2 = cat(1, curS(idCell2, 1:minNumTrial).matTS);
data1(:, iRun) = curTS1;
data2(:, iRun) = curTS2;
end

movingwin = [10 0.5];
params.taper = [10*6 119];
[C,phi,S12,S1,S2,t,f]=cohgramc(data1,data2,movingwin,params);
data2_s = data2(:, [2:end, 1]);
[sC,phi,S12,S1,S2,t,f]=cohgramc(data1,data2_s,movingwin,params);
figure
imagesc((C-sC)')