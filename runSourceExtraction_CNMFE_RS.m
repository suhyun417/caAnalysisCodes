% runSourceExtraction_CNMFE_RS.m
%
% Script;
% Run CNMF-E algorithm on resting state data from multiple sessions using "doCNMFE_1p.m"
% function (modified from demo_large_data_1p.m from CNMF-E toolbox)

clear all;

addpath('/projects/parksh/_toolbox/CNMF_E/');
cnmfe_setup;

% Set the directory
nameSubj = 'Tabla'; %'Max'; %'Tabla'; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; % 'Max'; % 'Tabla';
% setDateSession = {'20191113', '20191114', '20191118'}; %'20191223'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty

switch lower(nameSubj)
    case 'tabla'
        startSession = 1;
    case 'max'
        startSession = 1; %2;
end
 
for iSession = 1:5 %length(setDateSession) %startSession:length(setDateSession)

    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    % filename
    fname = fullfile(dirPreproc, 'ConcatRuns_RS.tif');
    [p, n, ext] = fileparts(fname);
    
    if isempty(dir(fullfile(p, [n '_source_extraction'], 'Sources2D_RS*.mat'))) %~exist(fullfile(p, [n '_source_extraction']), 'dir')
        
        paramCNMFE.gSig = 5; %3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
        paramCNMFE.gSiz = 13;          % pixel, neuron diameter
        
        fname_mat = fullfile(p, [n '.mat']); %'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/Session/20191119/_preproc/ConcatRuns_BPM_DFL.mat'
        matObj = matfile(fname_mat);
        Y = matObj.Y(1:matObj.Ysiz(1,1), 1:matObj.Ysiz(2,1), 1:1000);
        
        options.d1 = size(Y,1);
        options.d2 = size(Y,2);
        options.gSig = paramCNMFE.gSig;
        options.gSiz = paramCNMFE.gSiz;
        options.center_psf = 1;
        
        fprintf(1, 'Computing Cn and PNR for Session %d/%d: %s...\n', iSession, length(setDateSession), dateSession);
        [Cn, PNR] = correlation_image_endoscope(Y, options);
        sortPNR = sort(PNR(:));
        sortCn = sort(Cn(:));
        
        paramCNMFE.critPNR = 0.95;
        paramCNMFE.sortPNR = sortPNR;
        paramCNMFE.min_pnr = sortPNR(round(length(sortPNR).*paramCNMFE.critPNR));
        paramCNMFE.critCn = 0.95;
        paramCNMFE.sortCn = sortCn;
        paramCNMFE.min_corr = sortCn(round(length(sortCn).*paramCNMFE.critCn));
        
%         critPNR.crit95(iSession,1) = sortPNR(round(length(sortPNR).*0.95));
%         critPNR.crit90(iSession,1) = sortPNR(round(length(sortPNR).*0.90));
%         critCn.crit95(iSession,1) = sortCn(round(length(sortCn).*0.95));       
%         critCn.crit90(iSession,1) = sortCn(round(length(sortCn).*0.90));
        
        % parameters
        paramCNMFE.ssub = 1;
        paramCNMFE.memory_size_to_use = 30; %  64, ... % GB, memory space you allow to use in MATLAB
        paramCNMFE.memory_size_per_patch = 3; % ... 0.5, ...   % GB, space for loading data within one patch
        paramCNMFE.patch_dims = [64, 64]; %,...  %GB, patch size
%         paramCNMFE.batch_frames = 8000;
        paramCNMFE.Fs = 10;
        paramCNMFE.tsub = 1;
        paramCNMFE.K = [];
%         paramCNMFE.min_corr = 0.9; %0.85; %0.8; %0.9; %0.9;     % minimum local correlation for a seeding pixel
%         paramCNMFE.min_pnr = 20; %10;
        
        [neuron, flags] = doCNMFE_1p(fname, paramCNMFE);
        
        neuron.compress_results();
        paramCNMFE.flags = flags;
        
        
        [p, n, ext] = fileparts(fname);
        neuron.P.log_folder = fullfile(p, [n '_source_extraction']);
        file_path = fullfile(neuron.P.log_folder,  ['Sources2D_RS_', strrep(get_date(), ' ', '_'), '.mat']);
        
        save(file_path, 'neuron', 'paramCNMFE', '-v7.3');
        
        %     neuron.save_workspace_batch();
        
        clear neuron paramCNMFE mat_data*
        
    end
    
end