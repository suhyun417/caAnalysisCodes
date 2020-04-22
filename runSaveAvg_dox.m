% runSaveAvg_dox.m
%
% For each daily session and each resting state recording, compute average
% frame (across time) and save it as .mat and .tiff
% 2020/04/22 SHP


clear all;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/TIFFstack');
addpath('/projects/parksh/_toolbox/Fast_Tiff_Write/');
addpath('/projects/parksh/_toolbox/imagetools/');

%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

flagSaveFigure = 1;
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

for iSubj = 1:length(setSubj)
    nameSubj = setSubj{iSubj};
    
    % get session info
    [infoSession, opts] = readInfoSession_dox(nameSubj);

    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    clear infoSession
    
    for iSession = 1:nSession
        dateSession = setDateSession{iSession};

        dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        
        % get the info for concatenated runs & exp types
        load(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_all.mat'), 'paramConcat');       
        infoSession = paramConcat.infoSession;
        
        fprintf(1, ' %s: Session #%d/%d (%s): computing average frame run-by-run... \n', nameSubj, iSession, nSession, dateSession);
        
        for iRun = 1:length(infoSession)            
            
            load(fullfile(dirProcdata_session, '_preproc', sprintf('%s_sDS_cat.mat', infoSession(iRun).ImagingFilename)), 'Yf_cat')
            avgY = mean(Yf_cat, 3); % average across time dimension
            varY = var(Yf_cat, 0, 3); % compute variance across time dimension
            
            name_avg = sprintf('RS%d_Avg_%s_sDS_cat', iRun, infoSession(iRun).ImagingFilename);
            name_var = sprintf('RS%d_Var_%s_sDS_cat', iRun, infoSession(iRun).ImagingFilename);
            
            % save .mat file
            save(fullfile(dirProcdata_session, [name_avg, '.mat']), 'avgY');
            save(fullfile(dirProcdata_session, [name_var, '.mat']), 'varY');
            
            % save .tif file
            fastTiffStackWrite(fullfile(dirProcdata_session, [name_avg, '.tif']), avgY);
            fastTiffStackWrite(fullfile(dirProcdata_session, [name_var, '.tif']), varY);
                       
            fprintf(1, '               :: Run #%d/%d (%s): Average and variance across time are saved as .mat and .tif \n', ...
                iRun, length(infoSession) , infoSession(iRun).ImagingFilename);
        
        end
        
    end % session
    
end % subject


               
        
        
        
        