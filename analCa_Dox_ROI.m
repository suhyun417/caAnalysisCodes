% analCa_Dox_ROI.m
%
% analyze calcium data following dox administration based on ROIs
% using averaged frame

clear all;

%% Setting
% addpath('/projects/parksh/_toolbox/TIFFstack');
% addpath('/projects/parksh/_toolbox/Fast_Tiff_Write/');
% addpath('/projects/parksh/_toolbox/imagetools/');
addpath('/projects/parksh/_toolbox/ReadImageJROI/');


flagSaveFigure = 1;
dirFig = '/projects/ntis/figs/';

flagSaveFile = 1;


%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

for iSubj = 1:length(setSubj)
    nameSubj = setSubj{iSubj};
    
    % get session info
    [infoSession, opts] = readInfoSession_dox(nameSubj);

    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    clear infoSession
    
    % dox directory for each subject
    dirDox = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Dox', nameSubj);
    % Read ROI data
    d_roi = dir(fullfile(dirDox, '*.zip'));
    [setROIs] = ReadImageJROI(fullfile(d_roi.folder, d_roi.name));
    
    resultswholeFOV = struct([]);
    resultsROI = struct([]);
    for iSession = 1:nSession
        dateSession = setDateSession{iSession};

        dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        
        % retrieve file names of average frames
        d_avg = dir(fullfile(dirProcdata_session, '*Avg*.mat'));        
        
        meanF_FOV_run = NaN(length(d_avg), 1);
        meanF_ROI = NaN(length(d_avg), length(setROIs)); % run x roi
        FOV_ROI = cell(length(d_avg), length(setROIs));
        for iRun = 1:length(d_avg)            
            
            load(fullfile(d_avg(iRun).folder, d_avg(iRun).name))
            
            %% First, the entire FOV average
            meanF_FOV_run(iRun, 1) = mean(mean(avgY)); 
            
            for iRoi = 1:length(setROIs)
                % get the ROI coordinates
                coordsROI = setROIs{iRoi}.vnRectBounds; % [nTop nLeft nBottom nRight]
                
                % Compute the averaged fluorescence for each ROI
                FOV_ROI{iRun, iRoi} = avgY(coordsROI(1):coordsROI(3), coordsROI(2):coordsROI(4));
                meanF_ROI(iRun, iRoi) = mean(mean(FOV_ROI{iRun, iRoi}));
            end

        end % run
        
        resultswholeFOV(iSession).meanF_run = meanF_FOV_run; % vector of meanF over the entire FOV across runs
        resultswholeFOV(iSession).meanF = mean(meanF_FOV_run);
        
        resultsROI(iSession).setROIs = setROIs;
        resultsROI(iSession).infoRuns = d_avg;
        resultsROI(iSession).FOV_ROI_run = FOV_ROI;
        resultsROI(iSession).meanF_ROI_run = meanF_ROI;
        resultsROI(iSession).meanF_ROI = mean(meanF_ROI, 1);
        
        paramSession(iSession).dateSession = dateSession;
        paramSession(iSession).infoRuns = d_avg;
        
        
        if flagSaveFile
            fname_save = sprintf('doxResults_FOVROI_%s.mat', nameSubj);
            save(fullfile(dirDox, fname_save), 'results*', 'paramSession');
        end
        
    end
    
end
        
        
        
        