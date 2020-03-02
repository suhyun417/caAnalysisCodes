function [neuron, flags] = runCNMFE_batch_1p(fname, paramCNMFE)


%% clear the workspace and select data
% clear; clc; close all;
addpath('/projects/parksh/_toolbox/CNMF_E/');
cnmfe_setup;

%% choose multiple datasets or just one
neuron = Sources2D();
% nams = select_files();
nams = {fname}; %{fullfile(dirPreproc, 'ConcatRuns_BPM_DFL.tif')}; %{'./data_1p.tif'};          % you can put all file names into a cell array; when it's empty, manually select files
nams = neuron.select_multiple_files(nams);  %if nam is [], then select data interactively

%% parameters
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', paramCNMFE.memory_size_to_use, ... %30, ...  64, ... % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', paramCNMFE.memory_size_per_patch, ... %3, ... 0.5, ...   % GB, space for loading data within one patch
    'patch_dims', paramCNMFE.patch_dims, ... %[50, 50],...  %[64, 64],...  %GB, patch size
    'batch_frames', paramCNMFE.batch_frames);  %8000); %1000);           % number of frames per batch
% -------------------------      SPATIAL      -------------------------  %
gSig = paramCNMFE.gSig; %5; %3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = paramCNMFE.gSiz; %13;          % pixel, neuron diameter
ssub = paramCNMFE.ssub; %;           % spatial downsampling factor
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
Fs = paramCNMFE.Fs; %10; %15; %10;             % frame rate
tsub = paramCNMFE.tsub; %4; %2; %1; %5; %1;           % temporal downsampling factor
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
K = paramCNMFE.K; %[]; %10; %20; %[];             % maximum number of neurons per patch. when K=[], take as many as possible.
min_corr = paramCNMFE.min_corr; %0.9; %0.85; %0.8; %0.9; %0.9;     % minimum local correlation for a seeding pixel
min_pnr = paramCNMFE.min_pnr; %20; %10; %10; %10; %5; %8;       % minimum peak-to-noise ratio for a seeding pixel
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
min_corr_res = 0.85; %0.7;
min_pnr_res = 10; %6;
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
neuron.initComponents_batch(K, save_initialization, use_parallel);
% if show_init
%     figure();
%     ax_init= axes();
%     imagesc(Cn, [0, 1]); colormap gray;
%     hold on;
%     plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
% end

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

flags.save_demixed  =   save_demixed;
flags.save_initialization  = save_initialization;
flags.show_init = show_init;
flags.show_merge = show_merge;
flags.use_parallel = use_parallel;
flags.with_dendrites = with_dendrites;
flags.with_manual_intervention =with_manual_intervention;

% %% save workspace
neuron.save_workspace_batch();

