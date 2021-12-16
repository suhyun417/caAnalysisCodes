function [] = doAssembleCellCenterInfoAcrossSessions()

% 2021/12/14 SHP
% Aseemble cell center positions across daily sessions for a given FOV
% using the shifts saved in SUBJ_FOV#_shifts.mat
% - Create a grid space in n x n pixel resolution
% Pool in the pixel space. imagine a grid in n x n pixel
% resolution (e.g. 3 pixels? I can choose the criterion) and make a list of
% grid in columnar way. Then gather cell IDs from the REF and SES 

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
    
    
    
    
    
    
    
    

end % end of function