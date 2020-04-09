% runSourceExtraction_CNMFE_dox.m
%
% Script;
% Run CNMF-E algorithm for multiple sessions using a function called "runCNMFE_1p.m"
%   - modified from doSourceExtraction_CNMFE.m and originally from demo_batch_1p.m from CNMF-E toolbox
% 2020/04/09 SHP

clear all;

addpath('/projects/parksh/_toolbox/CNMF_E/');
cnmfe_setup;

% Set the directory
nameSubj = 'Tabla'; %'Max'; 

[infoSession, opts] = readInfoSession_dox(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty

clear infoSession

fid = fopen(fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', 'infoDoxSessions_minPNRCn.txt'), 'w');
for iSession = 1:length(setDateSession) %startSession:length(setDateSession)

    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    % filename
    fname = fullfile(dirPreproc, 'ConcatRuns_all.mat'); %fullfile(dirPreproc, 'ConcatRuns_BPM_DFL.tif');
    [p, n, ext] = fileparts(fname);
    
    if isempty(dir(fullfile(p, [n '_source_extraction'], 'Sources2D_all*.mat'))) %~exist(fullfile(p, [n '_source_extraction']), 'dir')
        
        paramCNMFE.gSig = 5; %3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
        paramCNMFE.gSiz = 13;          % pixel, neuron diameter
        
        fname_mat = fullfile(p, [n '.mat']); %'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/Session/20191119/_preproc/ConcatRuns_BPM_DFL.mat'
        matObj = matfile(fname_mat);
        setFrame = 1:1000;
        if strcmp(dateSession, '20191114')
            setFrame = 6001:7000;
        end
        Y = matObj.Y(1:matObj.Ysiz(1,1), 1:matObj.Ysiz(2,1), setFrame);
        
        options.d1 = size(Y,1);
        options.d2 = size(Y,2);
        options.gSig = paramCNMFE.gSig;
        options.gSiz = paramCNMFE.gSiz;
        options.center_psf = 1;
        
        fprintf(1, 'Computing Cn and PNR for Session %d/%d: %s...\n', iSession, length(setDateSession), dateSession);
        [Cn, PNR] = correlation_image_endoscope(Y, options);
        sortPNR = sort(PNR(:));
        sortCn = sort(Cn(:));
        
        paramCNMFE.critPNR = 0.9; %0.95;
        paramCNMFE.sortPNR = sortPNR;
        paramCNMFE.min_pnr = sortPNR(round(length(sortPNR).*paramCNMFE.critPNR));
        paramCNMFE.critCn = 0.9; %0.95;
        paramCNMFE.sortCn = sortCn;
        paramCNMFE.min_corr = sortCn(round(length(sortCn).*paramCNMFE.critCn));
        fprintf(1, 'Session %d/%d (%s): min_pnr = %2.2f, min_corr = %2.2f ...\n', ...
            iSession, length(setDateSession), dateSession, paramCNMFE.min_pnr, paramCNMFE.min_corr);
        fprintf(fid, 'Session %d/%d (%s): min_pnr = %2.2f, min_corr = %2.2f ...\n', ...
            iSession, length(setDateSession), dateSession, paramCNMFE.min_pnr, paramCNMFE.min_corr);
        
        clear Y Cn PNR sort*
        
%         critPNR.crit95(iSession,1) = sortPNR(round(length(sortPNR).*0.95));
%         critPNR.crit90(iSession,1) = sortPNR(round(length(sortPNR).*0.90));
%         critCn.crit95(iSession,1) = sortCn(round(length(sortCn).*0.95));       
%         critCn.crit90(iSession,1) = sortCn(round(length(sortCn).*0.90));
        
        % parameters
        paramCNMFE.ssub = 1;
        paramCNMFE.memory_size_to_use = 30; %  64, ... % GB, memory space you allow to use in MATLAB
        paramCNMFE.memory_size_per_patch = 3; % ... 0.5, ...   % GB, space for loading data within one patch
        paramCNMFE.patch_dims = [50, 50]; %,...  %[64, 64],...  %GB, patch size
%         paramCNMFE.batch_frames = 8000;
        paramCNMFE.Fs = 10;
        paramCNMFE.tsub = 1; %4;
        paramCNMFE.K = [];
%         paramCNMFE.min_corr = 0.9; %0.85; %0.8; %0.9; %0.9;     % minimum local correlation for a seeding pixel
%         paramCNMFE.min_pnr = 20; %10;
        
%         [neuron, flags] = runCNMFE_batch_1p(fname, paramCNMFE);
        [neuron, flags] = runCNMFE_1p(fname, paramCNMFE);
        
        neuron.compress_results();
        paramCNMFE.flags = flags;
        
        
        [p, n, ext] = fileparts(fname);
        neuron.P.log_folder = fullfile(p, [n '_source_extraction']);
        file_path = fullfile(neuron.P.log_folder,  ['Sources2D_all_nobatch_residualOn_customCritPNRCn0p9', strrep(get_date(), ' ', '_'), '.mat']);
        
        save(file_path, 'neuron', 'paramCNMFE', '-v7.3');
        
        %     neuron.save_workspace_batch();
        
        clear neuron paramCNMFE mat_data*
        
    end
    
end        
fclose(fid);

