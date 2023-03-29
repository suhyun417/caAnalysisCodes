% analCa_RS_PCA.m
%
% 2023/3/28 SHP modified from analCa_DFL_PCA.m
% perform PCA on resting state (spontaneous activity) responses and save the
% results of PCA

clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/PROJECTS/parksh';
        dirProcdata = '/Volumes/NIFVAULT/PROCDATA/parksh';
        dirRawdata = '/Volumes/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');


%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

for iSubj = 1:2 %2; %1;
    
    nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
    FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
    [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    
    %% load saved files
    % cells pooled across days
    fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat',...
        nameSubj, FOV_ID, nameSubj, FOV_ID));
    load(fname_stack, 'cellIDAcrossDay'); %, 'stackCellCenter')
    
    % cell quality info
    fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
        nameSubj, FOV_ID, nameSubj, FOV_ID));
    load(fname_cellQC, 'infoCells')
    
    % translational shift across days
    fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
        nameSubj, FOV_ID, nameSubj, FOV_ID));
    load(fname_shifts, 'shifts')
    
    % aligned cells TS 
    clear cellTS_RS resultsRS 
    fname_caTSFOV_RS = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_RSsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    load(fname_caTSFOV_RS, 'cellTS_RS', 'resultsRS')
    
    resultsPCA_session = struct([]);
    resultsPCA_2min = struct([]);
    matTS = struct([]);
    for iSession = 1:length(setDateSession)
        
        [coeffRun, scoreRun, latentRun, tsquaredRun, explainedRun] = pca(zscore(resultsRS(iSession).C_raw, 0, 2)); %pca(zscore(matAvgTS1)');
        [sortedScore, indCellRun] = sort(scoreRun(:,1:10), 'descend');
        
        resultsPCA_session(iSession).explained = explainedRun;
        resultsPCA_session(iSession).coeff = coeffRun(:, 1:10);
        resultsPCA_session(iSession).score = scoreRun(:, 1:10);
        resultsPCA_session(iSession).indCellSorted = indCellRun;
        
        curMatTS = resultsRS(iSession).C_raw;
        nSplit = size(curMatTS, 2)./1200;
        curMatTS_rs = reshape(curMatTS, size(curMatTS,1), size(curMatTS, 2)/nSplit, nSplit);
        curMatTS_rs_z = zscore(curMatTS_rs, 0, 2);
        
        for iSplit = 1:nSplit % 2 minute chunks
            
            [coeff, score, latent, tsquared, explained] = pca(squeeze(curMatTS_rs_z(:,:,iSplit))); %pca(zscore(matAvgTS1)');
            [sortedScore, indCell] = sort(score(:,1:10), 'descend');
                      
            resultsPCA_2min(iSession, iSplit).explained = explained;
            resultsPCA_2min(iSession, iSplit).coeff = coeff(:, 1:10);
            resultsPCA_2min(iSession, iSplit).score = score(:, 1:10);
            resultsPCA_2min(iSession, iSplit).indCellSorted = indCell;
            
%             figure;
%             imagesc(curMatTS_rs_z(resultsPCA_2min(iSession, iSplit).indCellSorted(:,1),:,iSplit), [-1 10])
            
        end
        
        matTS(iSession).matTS_session_raw = curMatTS;
        matTS(iSession).matTS_session_z = zscore(curMatTS, 0, 2);
        matTS(iSession).matTS_2min_raw = curMatTS_rs;
        matTS(iSession).matTS_2min_z = curMatTS_rs_z;
        
    end
    resultsPCA(iSubj).resultsPCA_session = resultsPCA_session;
    resultsPCA(iSubj).resultsPCA_2min = resultsPCA_2min;
    
    paramPCA(iSubj).matTS = matTS;
    
    
end


save(fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', 'RS_TS_PCA.mat'), 'resultsPCA', 'paramPCA')
    
    
    
    