% split the time series from the Sources2D_all*.mat file from CNMF-E output
% and save them to different experiment type

% getSessionInfo.m

clear all; close all;

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


%% directory
nameSubj = 'Tabla'; % 'Max'; % 'Tabla'; %'Max'; %'Tabla'; % 'Max'; % 'Tabla';
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

clear infoSession

for iSession = 2%1:length(setDateSession)

dateSession = setDateSession{iSession}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');


%% retrieve the experiment name for each run and save each source extracted time series
% get the info for concatenated runs & exp types
load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_all.mat'), 'paramConcat');
infoSession = paramConcat.infoSession;
% fname = sprintf('info_%s_%s.txt', dateSession, nameSubj); % example text file containing ename info for a given session
% fileID = fopen(fname);
% c = textscan(fileID, '%s %s %s', 'CommentStyle', '##'); 
% fclose(fileID);
    
% % get session info
% [infoSession, opts] = readInfoSession(nameSubj);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);

%% Load the source extraction data
addpath(fullfile(dirProjects, '_toolbox/CNMF_E/')); 
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
if isempty(d_sources2D)
    copyfile(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_all_source_extraction/Sources2D_all*.mat'), dirProcdata_session)
end
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); 


%% 
% order of concatenation used for source extraction

% listRun = cellstr(cat(1, infoSession.ImagingFilename));
% d_all = dir(fullfile(dirPreproc, '*_sDS_cat.mat'));
% listRun = {d_all.name}';
% % listRun_time = regexp(listRun, '\d{6}(?=_sDS)', 'match'); % Runs that are concatenated and used
% % listRun_time = vertcat(listRun_time{:});

count = 0;
for iRun = 1:length(infoSession)
    
%     nameRun = infoSession(iRun).ImagingFilename; %listRun{iRun}; %regexp(listRun{iRun}, '\d{6}(?=_sDS)', 'match'); %
%     ename = c{2}{contains(c{1}, nameRun)};
    
    % get the # of frames from each session without loading the matrix
    matObj = matfile(fullfile(dirPreproc, sprintf('%s_sDS_cat.mat', infoSession(iRun).ImagingFilename))); %listRun{iRun}));
    [d1 d2 nFrame] = size(matObj, 'Yf_cat');
    
%     infoRun(iRun) = infoSession(iRun);
    tSeries(iRun).idNeuron = neuron.ids;
    tSeries(iRun).C_raw = neuron.C_raw(:, count+1:count+nFrame);  %
    tSeries(iRun).C = neuron.C(:, count+1:count+nFrame);  %
    
    spatial(iRun).A = neuron.A;
    spatial(iRun).batches = neuron.batches;   
    
    count = count + nFrame;    
end

% save(fullfile(dirProcdata_session, 'all_ts.mat'), 'infoSession', 'tSeries', 'spatial')
clear catEName
 [catEName{1:length(infoSession), 1}] = deal(infoSession.ExpName);

locBPM = contains(catEName, 'BPM');
tSeries_BPM = tSeries(locBPM);
spatial_BPM = spatial(locBPM);
infoRun_BPM = infoSession(locBPM);
save(fullfile(dirProcdata_session, 'BPM_ts.mat'), 'infoRun_BPM', 'tSeries_BPM', 'spatial_BPM')

locDFL = contains(catEName, 'DFL');
tSeries_DFL = tSeries(locDFL);
spatial_DFL = spatial(locDFL);
infoRun_DFL = infoSession(locDFL);
save(fullfile(dirProcdata_session, 'DFL_ts.mat'), 'infoRun_DFL', 'tSeries_DFL', 'spatial_DFL')

locRS = contains(catEName, 'RS');
tSeries_RS = tSeries(locRS);
spatial_RS = spatial(locRS);
infoRun_RS = infoSession(locRS);
save(fullfile(dirProcdata_session, 'RS_ts.mat'), 'infoRun_RS', 'tSeries_RS', 'spatial_RS')

end

% % Using adjusted (ring background subtracted) time series from CNMFe
%     count = 0;
%     for iRun = 1:length(listRun_BPM)
%         
%         % get the # of frames from each session without loading the matrix
%         matObj = matfile(fullfile(dirProcdata_session, d_file_imaging(iRun).name));
%         [d1 d2 nFrame] = size(matObj, 'Y');
% 
%         for iCell = 1:size(neuron.A, 2)
% %             clear tempTS
% %             tempTS = matDFF(neuron.A(:,iCell)>0, :);
%             
%             tSeries_BPM(iCell, iRun).idNeuron_org =neuron.ids(iCell);
% %             tSeries_BPM(iCell, iRun).matTS = tempTS;
% %             tSeries_BPM(iCell, iRun).nPixel = size(tempTS, 1);
%             tSeries_BPM(iCell, iRun).mnTS_C_raw = neuron.C_raw(iCell, count+1:count+nFrame);  %mean(tempTS);
%             tSeries_BPM(iCell, iRun).mnTS_C = neuron.C(iCell, count+1:count+nFrame);
% %             tSeries_BPM(iCell, iRun).steTS = std(tempTS)./(sqrt(size(tempTS,1))-1);
%             
%         end   
%         count = count+nFrame;
%     end


    
    
% listRun_BPM = c{1}(contains(c{2}, 'BPM'));
% listRun_DFL = c{1}(contains(c{2}, 'DFL'));
% listRun_RS = c{1}(contains(c{2}, 'RS'));



% list of BPM imaging files
% listRun_BPM = {'103640', '104434', '105216'}'; % for 20191125_Tabla % eventually you want to retrieve this info from separate log file or something
% listRun_BPM = {'122312', '123102', '123609'}'; % for 20191125_Max
% listRun_BPM = {'123055', '123716', '124400', '125049', '125746', '130541'}'; % for 20191113_Tabla
% listRun_BPM = {'110653', '111411'}'; %{'105847', '110653', '111411'}'; %20191113_Max
