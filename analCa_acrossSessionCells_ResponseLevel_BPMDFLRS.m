% analCa_acrossSessionCells_ResponseLevel_BPMDFLRS.m
%
% 2024/06/04 SHP
% For each neuron, compute the response level (mean response before
% normalization) and noise level for each condition separately
%  -use the cell registration matrix to perform the anlaysis on the entire
%  population
%  -save the results as "SUBJ_FOV#_cellAcrossDay_responseLevel.mat"


clearvars;


%% Directory settings
directory = setDir_shp;
dirProjects = directory.dirProjects;
dirProcdata = directory.dirProcdata;
dirRawdata = directory.dirRawdata;
dirFig = directory.dirFig;


addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));

flagSaveFile = 0; %1; %0; %1;



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
    
    %% load DFL and BPM and RS in the workspace?
    resultsDFL = struct([]);
    for iS = 1:length(setDateSession)
        dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
        dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        load(fullfile(dirProcdata_session, 'DFL_ts_tML'));

        resultsDFL(iS).tS_session = tS_session;
    end
    % noise during the DFL: tS_session(1).matTS_C_raw-tS_session(1).matTS_C
    % % for all the runs

    resultsBPM = struct([]);
    for iS = 1:length(setDateSession)
        dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
        dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        load(fullfile(dirProcdata_session, 'BPM_ts_tML'));

        %     resultsBPM(iS).tS_session = tS_session;
        resultsBPM(iS).tS_session_stim = tS_session_stim;

        %     tS_session.tS_trial(iCell, iTrial)
        %     tS_sesion_stim(iCell, iStim) %% cells x 25 conditions
    end
    load(fullfile(dirProcdata_session, 'BPM_ts_tML.mat'), 'tS_run', 'tS_session') % what's the difference?
    % prepare the cell traces
    matTS_BPM = [];
    for iCell = 1:length(setCell)
        matTS_BPM(:, iCell) = cat(1, tS_run(1).tS_trial(setCell(iCell),:).matTS); % this is C_raw
    end
%     tS_session(1).tS_trial = cat(2, tS_run.tS_trial); % Cell by Trial



    % is noise level similar to all the conditions
    % first check the PNRs consistency across runs
    load(fullfile(dirProcdata_session, 'DFL_ts_tML.mat'))
    pnrs_dfl_set = squeeze(max(tS_session(1).matTS_C)./std(tS_session(1).matTS_C_raw-tS_session(1).matTS_C)); % for all the runs

%     load(fullfile(dirProcdata_session, 'BPM_ts.mat'))
%     for iRun = 1:length(tSeries_BPM)
%         pnrs_bpm_set(:,iRun) = max(tSeries_BPM(1).C, [], 2)./std(tSeries_BPM(1).C_raw-tSeries_BPM(1).C, 0, 2); % p
%     end



    iRun = 1; iMovie = 1;
    matTS_DFL = [];
    matTS_DFL = tS_session(iMovie).matTS_C_raw(:, setCell, iRun);

    load(fullfile(dirProcdata_session, 'RS_ts.mat'))
    matTS_RS = [];
    matTS_RS = tSeries_RS(1).C_raw(setCell, :)';


    % final data structure
    % respBPM.mean(iCell, iSession) = x % does the mean activity in this case make sense?
    % respBPM.max(iCell, iSession) = x
    % respBPM.var(iCell, iSession) = x
    % respBPM.ste(iCell, iSession) = x
    % respBPM.varNoise(iCell, iSession) = x
    % respBPM.nT(iCell, iSession) =length(xx); %
    
    % respDFL.max(iCell, iSession) = x
    % respDFL.var(iCell, iSession) = x
    % respDFL.ste(iCell, iSession) = x
    % respDFL.varNoise(iCell, iSession) = x
    % respDFL.nT(iCell, iSession) =length(xx); %

    % respRS.max(iCell, iSession) = x
    % respRS.var(iCell, iSession) = x
    % respRS.ste(iCell, iSession) = x
    % respRS.varNoise(iCell, iSession) = x
    % respRS.nT(iCell, iSession) =length(xx); %




    
    %%
    cellTS = struct([]);
    cellPix = struct([]);
    
    for iCell = 1:size(cellIDAcrossDay, 1)
        curCells_session = find(~isnan(cellIDAcrossDay(iCell, :)));
        curCells_id = cellIDAcrossDay(iCell, curCells_session);

        % for each valid session
        for iS = 1





    end


    %%
    if flagSaveFile
        fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
        save(fname_caTSFOV, 'cellTS', 'cellPix', 'resultsDFL')
        fprintf(1, '\n Saving files for %s FOV %d\n', nameSubj, FOV_ID)
    end

end

% snrs = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2); % signal variance divided by noise variance
% pnrs = max(neuron.C, [], 2)./std(neuron.C_raw-neuron.C, 0, 2); % peak amplitude divided by noise std
% varNoise = var(neuron.C_raw-neuron.C, 0, 2); % variance of noise