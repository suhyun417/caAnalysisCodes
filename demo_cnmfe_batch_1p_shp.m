%% clear the workspace and select data 
clear; clc; close all; 
addpath('/projects/parksh/_toolbox/CNMF_E/');
cnmfe_setup

%% Set the directory
nameSubj = 'Max'; % 'Tabla'; %'Max'; %'Tabla'; % 'Max'; % 'Tabla';
dateSession = '20191223'; %'20191113'; %'20191125';
% datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

dirProcdata_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession);
dirPreproc = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession, '_preproc');

% dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
% list of BPM ML files and imaging files  %: eventually you want to retrieve this info from separate log file instead of hard coding
% listRun_BPM = {'103640', '104434', '105216'}'; % for 20191125_Tabla 
% listRun_BPM = {'122312', '123102', '123609'}'; % for 20191125_Max
% listRun_BPM = {'123055', '123716', '124400', '125049', '125746', '130541'}'; %20191113_Tabla
% listRun_BPM = {'110653', '111411'}'; %{'105847', '110653', '111411'}'; %20191113_Max


% %% concatenate the BPM runs into a single matrix for CNMFe source extraction
% % and save it in *.mat file
% % d = dir(fullfile(dirPreproc, 'DFL_*.mat')); 
% d = dir(fullfile(dirPreproc, 'BPM_*.mat'));
% if isempty(d)
%     fprintf(1, 'Concatenated .mat file for BPM does not exist. Creating one now...\n')
%     
%     Y = [];
%     for iRun = 1:length(listRun_BPM)
%         fprintf(1, '      Loading BPM run #%d/%d: %s_sDS_cat_RigidMC.mat \n', iRun, length(listRun_BPM), listRun_BPM{iRun})
%         load(fullfile(dirPreproc, sprintf('%s_sDS_cat_RigidMC.mat', listRun_BPM{iRun})), 'Mr')
%         
%         % retrieve and concatenate in time
%         Y = cat(3, Y, Mr);
%         clear Mr
%     end
%     Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
%     
%     fprintf(1, 'Concatenation is done! The file is being saved as BPM_sDS_cat_RigidMC_cat.mat...\n')
%     save(fullfile(dirPreproc, 'BPM_sDS_cat_RigidMC_cat.mat'), 'Y', 'Ysiz', '-v7.3')
%     fprintf(1, '...Done!\n')
% end

% %% concatenate DFL runs into a single matrix for CNMFe source extraction
% % and save it in *.mat file
% listRun_DFL = {'131229', '131728', '132314', '132823', '133237', '133624'}; %20191113_Tabla
% 
% fname_catmat = 'DFL_sDS_cat_RigidMC_cat.mat';
% d = dir(fullfile(dirPreproc, fname_catmat)); 
% % d = dir(fullfile(dirPreproc, 'BPM_*.mat'));
% if isempty(d)
%     fprintf(1, 'Concatenated .mat file for DFL does not exist. Creating one now...\n')
%     
% %     Y = [];
% %     data = matfile(fullfile(dirPreproc, fname_catmat), 'Writable', true); 
%     count = 0; % for concatenation
%     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)    
%     for iRun = 1:length(listRun_DFL)
%         tic;
%         fprintf(1, '      Loading DFL run #%d/%d: %s_sDS_cat_RigidMC.mat \n', iRun, length(listRun_DFL), listRun_DFL{iRun})
%         load(fullfile(dirPreproc, sprintf('%s_sDS_cat_RigidMC.mat', listRun_DFL{iRun})), 'Mr')
%         
%         % write the file into a single mat file
%         [d1, d2, T] = size(Mr);
%         if iRun == 1
%             Y = Mr;
%             save(fullfile(dirPreproc, fname_catmat), 'Y', '-v7.3');
%         else 
%             data = matfile(fullfile(dirPreproc, fname_catmat), 'Writable', true); 
%             data.Y(:, :, count+1:count+T) = Mr;
%         end
% %         Y = cat(3, Y, Mr);
%         clear Mr        
%         count = count+T;
%         toc
%     end
%     Ysiz  = size(data, 'Y')'; % [d1, d2, T]'; % following CNMFe's mat file convention
%     save(fullfile(dirPreproc, fname_catmat), 'Ysiz', '-append'); %, '-v7.3')
% %     Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
%     
% %     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)
% %     save(fullfile(dirPreproc, fname_catmat), 'Y', 'Ysiz', '-v7.3')
%     fprintf(1, '...Done!\n')
% end
% 
% % fprintf('CNMF_E is converting TIFF file to *.mat file'); 
% % % create a mat file 
% % Tchunk = min(T, round(2^29/d1/d2)); %each chunk uses at most 4GB
% % Y = smod_bigread2(nam, 1, Tchunk);  %#ok<*NASGU>
% % save(nam_mat, 'Y', 'Ysiz', '-v7.3'); 
% % if Tchunk==T
% %     return; 
% % else
% %     data = matfile(nam_mat, 'Writable', true); 
% %     t0 = Tchunk+1; 
% %     while t0<=T
% %         num2read = min(t0+Tchunk-1, T) - t0 + 1; 
% %         tmpY = smod_bigread2(nam, t0, num2read); 
% %         data.Y(:, :, (1:num2read)+t0-1) = tmpY; 
% %         t0 = t0 + num2read; 
% %     end 
% % end 


% %% concatenate all runs into a single matrix for CNMFe source extraction
% % and save it in *.mat file
% % listRun_DFL = {'131229', '131728', '132314', '132823', '133237', '133624'}; %20191113_Tabla
% 
% fname_catmat = 'allRuns_sDS_cat_RigidMC_cat.mat';
% d = dir(fullfile(dirPreproc, fname_catmat)); 
% % d = dir(fullfile(dirPreproc, 'BPM_*.mat'));
% if isempty(d)
%     fprintf(1, 'Concatenated .mat file for all runs does not exist. Creating one now...\n')
%     
% %     Y = [];
% %     data = matfile(fullfile(dirPreproc, fname_catmat), 'Writable', true); 
%     count = 0; % for concatenation
%     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)   
%     
%     fname = sprintf('info_%s_%s.txt', dateSession, nameSubj); % example text file containing ename info for a given session
%     fileID = fopen(fname);
%     c = textscan(fileID, '%s %s %s', 'CommentStyle', '##');
%     fclose(fileID);
%     
%      listRun_cat = c{1}(contains(c{2}, {'BPM', 'DFL'}));
% %      d_all = dir(fullfile(dirPreproc, '*_sDS_cat_RigidMC.mat'));
%     
%      for iRun = 1:length(listRun_cat)
%         tic;
% %         fprintf(1, '      Loading run #%d/%d: %s \n', iRun, length(d_all), d_all(iRun).name)
% %         load(fullfile(dirPreproc, d_all(iRun).name), 'Mr'); %sprintf('%s_sDS_cat_RigidMC.mat', listRun_DFL{iRun})), 'Mr')
%         
%         fprintf(1, '      Loading run #%d/%d: %s \n', iRun, length(listRun_cat), listRun_cat{iRun})
%         load(fullfile(dirPreproc, sprintf('%s_sDS_cat_RigidMC.mat', listRun_cat{iRun})), 'Mr');
%         
%         % write the file into a single mat file
%         [d1, d2, T] = size(Mr);
%         if iRun == 1
%             Y = Mr;
%             save(fullfile(dirPreproc, fname_catmat), 'Y', '-v7.3');
%         else 
%             data = matfile(fullfile(dirPreproc, fname_catmat), 'Writable', true); 
%             data.Y(:, :, count+1:count+T) = Mr;
%         end
% %         Y = cat(3, Y, Mr);
%         clear Y Mr        
%         count = count+T;
%         toc
%     end
%     Ysiz  = size(data, 'Y')'; % [d1, d2, T]'; % following CNMFe's mat file convention
%     save(fullfile(dirPreproc, fname_catmat), 'Ysiz', '-append'); %, '-v7.3')
% %     Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
%     
% %     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)
% %     save(fullfile(dirPreproc, fname_catmat), 'Y', 'Ysiz', '-v7.3')
%     fprintf(1, '...Done!\n')
% end

%% choose multiple datasets or just one  
neuron = Sources2D(); 
% nams = select_files();
nams = {fullfile(dirPreproc, 'allRuns_sDS_cat_cat_RigidMC.tif')}; % {fullfile(dirPreproc, 'allRuns_sDS_cat_RigidMC_cat.mat')}; %'DFL_sDS_cat_RigidMC_cat.mat')}; %'BPM_sDS_cat_RigidMC_cat.mat')};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/20191125/_preproc/103640_sDS_cat_RigidMC.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Max/20191112/_preproc/tempCatAcTrial_MovieSet1_1.tif'}; %tempAvgAcTrial_MovieSet1_1.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Max/20191112/_preproc/112414_sDS_cat_RigidMC.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Max/20191112/_preproc/105441_sDS_cat_RigidMC.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/20191113/_preproc/131728_sDS_cat_RigidMC.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/20191113/_preproc/tempAvgAcTrial_MovieSet1_1.tif'};
% nams = {'/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/20191113/_preproc/123055_sDS_cat_RigidMC.tif'}; %{'./data_1p.tif'};          % you can put all file names into a cell array; when it's empty, manually select files 
nams = neuron.select_multiple_files(nams);  %if nam is [], then select data interactively 

%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 30, ...  64, ... % GB, memory space you allow to use in MATLAB 
    'memory_size_per_patch', 3, ... 0.5, ...   % GB, space for loading data within one patch 
    'patch_dims', [50, 50],...  %[64, 64],...  %GB, patch size 
    'batch_frames', 8000); %1000);           % number of frames per batch 
  % -------------------------      SPATIAL      -------------------------  %
gSig = 5; %3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = 13; %10; %13;          % pixel, neuron diameter
ssub = 1;           % spatial downsampling factor
with_dendrites = false; %true;   % with dendrites or not
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals';

% -------------------------      TEMPORAL     -------------------------  %
Fs = 10; %15; %10;             % frame rate
tsub = 4; %2; %1; %5; %1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

nk = 3;             % detrending the slow fluctuation. usually 1 is fine (no detrending) % 5 is following demo_endoscope.m
% when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 1.5; % 1.4;
ring_radius = round(bg_neuron_factor * gSiz);  % when the ring model used, it is the radius of the ring used in the background model.
%otherwise, it's just the width of the overlapping area
num_neighbors = 50; % number of neighbors for each neuron

% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step
merge_thr = 0.85; %0.65;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
merge_thr_spatial = [0.8, 0.4, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

% -------------------------  INITIALIZATION   -------------------------  %
K = []; %10; %20; %[];             % maximum number of neurons per patch. when K=[], take as many as possible.
min_corr = 0.9; %0.8; %0.9; %0.9;     % minimum local correlation for a seeding pixel
min_pnr = 15; %10; %10; %5; %8;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 1; %0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;    % use parallel computation for parallel computing
show_init = true;   % show initialization results
choose_params = true; % manually choose parameters
center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
% set the value as false when the background fluctuation is small (2p)

% -------------------------  Residual   -------------------------  %
min_corr_res = 0.7;
min_pnr_res = 6;
seed_method_res = 'auto';  % method for initializing neurons from the residual
update_sn = true;

% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false;

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not
kt = 3;                 % frame intervals

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', tsub, ...                       % -------- temporal --------
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization -----
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs;

%% distribute data and be ready to run source extraction 
neuron.getReady_batch(pars_envs); 

%% initialize neurons in batch mode 
[center, Cn, PNR] = neuron.initComponents_batch(K, save_initialization, use_parallel); 
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end

%% update background (%% originally after spatial and temporal, but here I'm following the demo_large_data_1p)
neuron.update_background_batch(use_parallel); 


%% udpate spatial components for all batches
neuron.update_spatial_batch_shp(use_parallel); %neuron.update_spatial_batch(use_parallel); 

%% udpate temporal components for all bataches
neuron.update_temporal_batch(use_parallel); 

%% update background 
neuron.update_background_batch(use_parallel); 

%% delete neurons 

%% merge neurons 

%% get the correlation image and PNR image for all neurons 
neuron.correlation_pnr_batch(); 

%% concatenate temporal components 
neuron.concatenate_temporal_batch(); 
% neuron.viewNeurons([],neuron.C_raw); 

%% save workspace 
neuron.save_workspace_batch(); 