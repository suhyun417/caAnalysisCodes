% split the time series from the Sources2D_all*.mat file from CNMF-E output
% and save them to different experiment type

% getSessionInfo.m

clear all; close all;

%% directory
nameSubj = 'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';


%% retrieve the experiment name for each run and save each source extracted time series
% get the info for concatenated runs & exp types
load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL.mat'), 'paramConcat');
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
addpath('/projects/parksh/_toolbox/CNMF_E/'); 
cnmfe_setup;
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
if isempty(d_sources2D)
    copyfile(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL_source_extraction/Sources2D_all*.mat'), dirProcdata_session)
end
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
load(fullfile(d_sources2D(1).folder, d_sources2D(1).name)); 


%% 
% order of concatenation used for source extraction
listRun = cellstr(cat(1, infoSession.ImagingFilename));
% d_all = dir(fullfile(dirPreproc, '*_sDS_cat.mat'));
% listRun = {d_all.name}';
% % listRun_time = regexp(listRun, '\d{6}(?=_sDS)', 'match'); % Runs that are concatenated and used
% % listRun_time = vertcat(listRun_time{:});

count = 0;
for iRun = 1:length(listRun)
    
    nameRun = regexp(listRun{iRun}, '\d{6}(?=_sDS)', 'match'); %
    ename = c{2}{contains(c{1}, nameRun)};
    
    % get the # of frames from each session without loading the matrix
    matObj = matfile(fullfile(dirPreproc, listRun{iRun}));
    [d1 d2 nFrame] = size(matObj, 'Mr');
    
    infoRun(iRun).nameRun = listRun{iRun};
    infoRun(iRun).ename = ename;
    tSeries(iRun).idNeuron = neuron.ids;
    tSeries(iRun).C_raw = neuron.C_raw(:, count+1:count+nFrame);  %
    tSeries(iRun).C = neuron.C(:, count+1:count+nFrame);  %
    
    count = count + nFrame;
    
end

save(fullfile(dirProcdata_session, 'all_ts.mat'), 'infoRun', 'tSeries')

catEName = cellstr(cat(1, infoRun.ename));
locBPM = contains(catEName, 'BPM');
tSeries_BPM = tSeries(locBPM);
infoRun_BPM = infoRun(locBPM);
save(fullfile(dirProcdata_session, 'BPM_ts.mat'), 'infoRun_BPM', 'tSeries_BPM')

locDFL = contains(catEName, 'DFL');
tSeries_DFL = tSeries(locDFL);
infoRun_DFL = infoRun(locDFL);
save(fullfile(dirProcdata_session, 'DFL_ts.mat'), 'infoRun_DFL', 'tSeries_DFL')

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
