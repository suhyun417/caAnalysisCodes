% runConcatAllRuns_dox.m
% 
% This script concatenates all the preprocessed files (*_mc.tif) for a
% given daily session and saves it as "ConcatRuns_all.mat" for further
% source extraction, using "doConcatRuns" function
%   - modified from "runConcatAllRuns.m"
% 2020/04/09 SHP

clear all;

%% for concatenation
% get session info
setNameSubj = {'Tabla', 'Max'};

for iSubj = 1:length(setNameSubj)

    nameSubj = setNameSubj{iSubj}; %'Tabla';
    [infoSession, opts] = readInfoSession_dox(nameSubj);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    for iSession = 1:nSession
        
        dateSession = setDateSession{iSession}; %'20191113'; 
        
        dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirProcdata_session_preproc = fullfile(dirProcdata_session, '_preproc');
        
        % List of runs to process: get the info from the xls file
        locRun_session = find((contains(infoSession.(1), dateSession).*(infoSession.(6)>0)));
        
        listRun = infoSession.(2)(locRun_session);
        listFileName = strcat(fullfile(dirProcdata_session_preproc, listRun), '*_mc.tif');
        % for iFile = 1:length(tempListFile)
        %     d = dir(tempListFile{iFile});
        %     listFileName{iFile, 1} = fullfile(;
        % end
        
        fname_cat = fullfile(dirProcdata_session_preproc, 'ConcatRuns_all');
        fprintf(1, 'Session %d/%d: %s: Concatenating all runs..\n', iSession, nSession, dateSession)
        tic; doConcatRuns(listFileName, fname_cat); toc;
        fprintf(1, 'Session %d/%d: %s: ............................Done!\n', iSession, nSession, dateSession)
        
        paramConcat.listRun = listRun;
        paramConcat.listFileName = listFileName;
        paramConcat.infoSession = table2struct(infoSession(locRun_session, :));
        
        save([fname_cat '.mat'], 'paramConcat', '-append')
        
        clear list* loc*
        
    end
end