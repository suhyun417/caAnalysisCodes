nameSubj = 'Max'; %'Tabla'; % 'Max'; % 'Tabla';
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

for iSession = 1:length(setDateSession)

dateSession = setDateSession{iSession}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
dirPreproc = fullfile(dirProcdata_session, '_preproc');

% delete a copied Sources2D_all*.mat in main "session" directory
d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
if ~isempty(d_sources2D)
    delete(fullfile(d_sources2D.folder, d_sources2D.name));
end

delete(fullfile(dirProcdata_session, '*ts.mat'))
delete(fullfile(dirProcdata_session, '*ts_tML.mat'))

% % delete the batch directories in "/_preproc/ConcatRuns_BPM_DFL_source_extraction"
% dirPreproc_SE = fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL_source_extraction');
% rmdir(fullfile(dirPreproc_SE, 'frames*'), 's')

% % make a new directory in "/_preproc/ConcatRuns_BPM_DFL_source_extraction" and move the previous Sources2D attempt to that directory
% movefile(fullfile(dirPreproc_SE, 'Sources2D_all*.mat'), fullfile(dirPreproc_SE, '_critCn0p9_critPNR0p9'))
% % fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL_source_extraction/Sources2D_all*.mat'), dirProcdata_session)

% % TEMPORAY SOLUTION: move the previous Sources2D_all file to main "session"
% % directory to try out other stuff while the source extraction happening
% copyfile(fullfile(dirPreproc_SE, '_minCorr0p9_minPnr20', 'Sources2D_all*.mat'), dirProcdata_session)

end