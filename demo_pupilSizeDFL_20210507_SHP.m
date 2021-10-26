% demo_pupilSizeDFL_20210507_SHP.m
%
% Read pupil size data from DFL session and start from there
%

setSubjName = {'Tabla', 'Max'};
iSubj = 1; %
% for iSubj = 1:length(setSubjName)
    clear Results*
    
    nameSubj =  setSubjName{iSubj}; %'Max'; %'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
    % dateSession = '20191113'; %'20191125';
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj);
    S = table2struct(infoSession);
    
    % setExpName = {S.ExpName}';
    setMLFilename = {S.MLFilename}';
    
    indDFLRuns = contains(setMLFilename, 'DFL') & cat(1, S.flagPreproc)>0; %% containing "DFL" in filename AND flagPreproc value of 1
    setFilename = setMLFilename(indDFLRuns);
    
    for iFile = 1:length(setFilename)
        filename = strcat(setFilename{iFile}, '.bhv2');
        
        dateSession = filename(1:6);
        
        if str2num(dateSession) < 191121
            dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
        else
            dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
        end
        
        % filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have
        
        %% Read the file
        data = mlread(fullfile(dirBHV, filename)); % mlread(filename);
        
        
    end
    
% end