function [shifts] = doRegisterSession(nameSubj, FOV_ID, flagSaveFile)
%
% 2021/12/13 SHP
% Performs longitudinal registration for a given subject and FOV 
% - Load reference image
% - For each following sessions
%   - load session image
%   - perform the rigid-body registration
%   - save the shifts (and registration parameters and other outcomes)
%       - for each animal and each FOV (e.g. Tabla_FOV1_shifts.mat)
% - INPUT
%   - nameSubj: full name of subject (string)
%   - FOV_ID: ID of FOV, in number (double)
%   - flagSaveFile: if 1, it will save the shifts in a .mat file
% - OUTPUT
%   - shifts: N by 2 array contains x, y shifts that need to be applied to
%   each of the N sessions



%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/projects/parksh';
        dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
        dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end
% dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

%% Session info & optional parameters

% get session info
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% Load the reference image (first session)
dirRefImage = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
    nameSubj, setDateSession{1}));
imgRef = loadtiff(fullfile(dirRefImage, 'mc_template.tif'));

shifts = NaN(nSession, 2);
shifts(1, :) = [0 0]; % 1st is always the ref
for iSession = 2:nSession
    % Load session image to align to the reference image
    dateSession = setDateSession{iSession};
    dirSessionImage = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
        nameSubj, dateSession));
    imgSession = loadtiff(fullfile(dirSessionImage, 'mc_template.tif'));
    
    %     paramHPF.gSig = 7;
    %     paramHPF.gSiz = 17;
    paramRegister.bound = 0; %40;
    %     paramHPF.imfilter = 'imfilter(Yf,psf,''symmetric'')';
    
    %     Yf = loadtiff(fname);
    [d1,d2,T] = size(imgSession);
    
    bound = paramRegister.bound; %40; %0;
    
    % set options
    %     [p, nameIn, ext] = fileparts(fname);
    clear options_r
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'max_shift',20,'iter',1,'correct_bidir',false, ...
        'output_type', 'mat'); %'tif', 'tiff_filename', './registrationTest_rgm.tif');%,'bin_width',T, ...
    paramRegister.options_r = options_r;
    
    % register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre(imgSession(bound/2+1:end-bound/2,bound/2+1:end-bound/2),options_r, imgRef); toc % register filtered data
    %      [M_final,shifts,template,options,col_shift] = normcorre(Y,options,template);
    % apply shifts and save it as desired format described in the options
    %     tic; Mrg = apply_shifts(imgSession,shifts1,options_r,bound/2,bound/2); toc
    
    shifts(iSession, :) = squeeze(shifts1.shifts)';
    diff(iSession, 1) = shifts1.diff;
    
end

if flagSaveFile
    fname_shifts = fullfile(dirProjects, sprintf('0Marmoset/Ca/tempData/%s_FOV%d_shifts.mat', nameSubj, FOV_ID));
%     fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    save(fname_shifts, 'shifts', 'diff', 'paramRegister')
end

% % Check the registration
% figure;
% subplot(1, 2, 1);
% imshowpair(imgRef, imgSession);
% title('Before Registration')
% subplot(1, 2, 2);
% imshowpair(imgRef, Mrg);
% title('Rigid Motion Regiration')
% % subplot(1, 3, 3);
% % imshowpair(imgRef, Mnrg);
% % title('Non-rigid motion Registration')

%         %% Apply shifts directly to neuron.A: yes you can just add the shifts
%         addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
%         cnmfe_setup;
%         
% %         % REF IMAGE
% %         d_sources2D_ref = dir(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
% %             nameSubj, setDateSession{1})));
% %         load(fullfile(d_sources2D_ref(1).folder, d_sources2D_ref(1).name));  
% % %         tName = ls(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
% % %             nameSubj, setDateSession{1}))); % because "dir" doesn't work 
% % %         load(strtrim(tName));
% %         imgFOV_ref = neuron.Cn.*neuron.PNR;
% %         [center_ref] = neuron.estCenter();
%                 
%         % SESSION IMAGE
%         clear neuron
%         dirProcdata_session = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/',...
%             nameSubj, dateSession));
%         d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));        
%         load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
% %         tName_ses = ls(fullfile(dirProcdata_session, 'Sources2D_all*'));
% %         load(strtrim(tName_ses))
%         imgFOV_ses = neuron.Cn.*neuron.PNR;
%         imgFOV_ses_reg = apply_shifts(imgFOV_ses, shifts1,options_r,bound/2,bound/2);
%         
% %         figure
% %         subplot(1,2,1);
% %         imshowpair(imgFOV_ref, imgFOV_ses)
% %         title('Before Registration')
% %         subplot(1,2,2);
% %         imshowpair(imgFOV_ref, imgFOV_ses_reg)
% %         title('After Registration')

%         [center_ses] = neuron.estCenter();
%         shifts = squeeze(shifts1.shifts);
%         center_ses_reg = center_ses + shifts';
        

%         figure
%         sp(1) = subplot(1,2,1);
%         plot(center_ref(:,2), center_ref(:, 1), 'r.')
%         hold on
%         plot(center_ses(:,2), center_ses(:, 1), 'b.')
%         sp(2) = subplot(1,2,2);
%         plot(center_ref(:,2), center_ref(:, 1), 'r.')
%         hold on
%         plot(center_ses_reg(:,2), center_ses_reg(:, 1), 'g.')
%         set(sp, 'YDir', 'reverse')
%         [d1,d2] = size(neuron.Cn);
%         set(sp, 'XLim', [0 d2], 'YLim', [0 d1])
%         title(sp(1), 'Cell center: before registration')
%         title(sp(2), 'Cell center: after registration')

end
