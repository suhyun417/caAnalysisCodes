% runAnalyzeML.m
%
% Script;
% For each imaging file, save timing information from corresponding bhv2 file
% Naming scheme: do we want to keep each imaging run name? or 

nameSubj =  'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
% dateSession = '20191113'; %'20191125';  
    
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

for iSession = 1:nSession
    
    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    if str2num(dateSession) < 20191121
        dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
    else
        dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
    end
    
    % get the info for concatenated runs & exp types
    load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_BPM_DFL.mat'), 'paramConcat');
    
    infoSession = paramConcat.infoSession;
    
    % process it run-by-run
    count_BPM = 0;
    count_DFL = 0;
    for iRun = 1:length(infoSession)
        
        fname_bhv = fullfile(dirBHV, [infoSession(iRun).MLFilename '.bhv2']);
        
        switch upper(infoSession(iRun).ExpName)
            case 'BPM'
                count_BPM = count_BPM+1;
                fname_mat = fullfile(dirProcdata_session, sprintf('BPM_%d_tML.mat', count_BPM));
                saveTimingML_BPM(fname_bhv, fname_mat);
                
                % additionally save relevant info with regard to imaging file
                infoImaging = infoSession(iRun);
                save(fname_mat, 'infoImaging', '-append')
            case 'DFL'
                count_DFL = count_DFL+1;
                fname_mat = fullfile(dirProcdata_session, sprintf('DFL_%d_tML.mat', count_DFL));
                saveTimingML_DFL(fname_bhv, fname_mat);
                
                % additionally save relevant info with regard to imaging file
                infoImaging = infoSession(iRun);
                save(fname_mat, 'infoImaging', '-append')
        end
    end
    
end

