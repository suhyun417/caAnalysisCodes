function [] = doAssembleCellCenterInfoAcrossSessions(nameSubj, FOV_ID, flagSaveFile)

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
addpath(fullfile(dirProjects, '_toolbox/CNMF_E/'));
cnmfe_setup;
% gcp; % for parallel processingls

%% Session info & optional parameters
% nameSubj = 'Tabla';
% FOV_ID = 1;
% get session info
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%%
fname_shifts = fullfile(dirProjects, sprintf('0Marmoset/Ca/tempData/%s_FOV%d_shifts.mat', nameSubj, FOV_ID));
load(fname_shifts, 'shifts')

stackCellCenter = [];
for iSession = 1:nSession

    d_sources2D = dir(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
        nameSubj, setDateSession{iSession})));
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
%     dims = size(neuron.Cn);
    [center] = neuron.estCenter() + shifts(iSession, :); % this is y and x for each cell in image coordinate ([0,0] is upper left corner)
    
    indValidCell = 1:length(center);
    if sum(max(round(center)) > size(neuron.Cn))|| sum(min(round(center)) <= 0)
        indValidCell = find(sum(cat(2, round(center(:,1))<size(neuron.Cn,1), round(center(:,2))<size(neuron.Cn,2), ...
            round(center(:,1)) > 0, round(center(:,2)) > 0 ), 2) > 3);
    end        
    
    locCell = sub2ind(size(neuron.Cn), round(center(indValidCell,1)), round(center(indValidCell,2)));
    curCanvas = NaN(size(neuron.Cn));
    curCanvas(locCell) = indValidCell;
    
    % figure
    % imagesc(curCanvas)
    % imagesc(curCanvas, 'AlphaData', ~isnan(curCanvas))
    % colormap(lines)
    
    stackCellCenter = cat(3, stackCellCenter, curCanvas);
end

% tempStack = sum(~isnan(stackCellCenter), 3); % # of cells per each pixel space
% figure;
% imagesc(tempStack) 
% 
% cat(2, squeeze(stackCellCenter(45, 219, :)), squeeze(stackCellCenter(42, 218, :)), squeeze(stackCellCenter(46, 218, :)), ...
%     squeeze(stackCellCenter(47, 218, :)), squeeze(stackCellCenter(44, 217, :)), squeeze(stackCellCenter(45, 217, :)))
% cat(2, squeeze(stackCellCenter(10, 105, :)), squeeze(stackCellCenter(10, 103, :)), squeeze(stackCellCenter(11, 104, :)), squeeze(stackCellCenter(9, 106, :)), squeeze(stackCellCenter(13, 104, :)))

if flagSaveFile
    fname_stack = fullfile(dirProjects, sprintf('0Marmoset/Ca/tempData/%s_FOV%d_stackedCenter.mat', nameSubj, FOV_ID));
%     fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_stackedCenter.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    save(fname_stack, 'stackCellCenter')
end



% critPixel = 2; % criterion to be near enough
% 
% 
% 
% 
% 
% shifts = NaN(nSession, 2);
% shifts(1, :) = [0 0]; % 1st is always the ref
% for iSession = 2:nSession
%     % Load session image to align to the reference image
%     dateSession = setDateSession{iSession};
%     dirSessionImage = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
%         nameSubj, dateSession));
%     imgSession = loadtiff(fullfile(dirSessionImage, 'mc_template.tif'));
%     
%     %     paramHPF.gSig = 7;
%     %     paramHPF.gSiz = 17;
%     paramRegister.bound = 0; %40;
%     %     paramHPF.imfilter = 'imfilter(Yf,psf,''symmetric'')';
%     
%     %     Yf = loadtiff(fname);
%     [d1,d2,T] = size(imgSession);
%     
%     bound = paramRegister.bound; %40; %0;
%     
%     % set options
%     %     [p, nameIn, ext] = fileparts(fname);
%     clear options_r
%     options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'max_shift',20,'iter',1,'correct_bidir',false, ...
%         'output_type', 'mat'); %'tif', 'tiff_filename', './registrationTest_rgm.tif');%,'bin_width',T, ...
%     paramRegister.options_r = options_r;
%     
%     % register using the high pass filtered data and apply shifts to original data
%     tic; [M1,shifts1,template1] = normcorre(imgSession(bound/2+1:end-bound/2,bound/2+1:end-bound/2),options_r, imgRef); toc % register filtered data
%     %      [M_final,shifts,template,options,col_shift] = normcorre(Y,options,template);
%     % apply shifts and save it as desired format described in the options
%     %     tic; Mrg = apply_shifts(imgSession,shifts1,options_r,bound/2,bound/2); toc
%     
%     shifts(iSession, :) = squeeze(shifts1.shifts)';
%     diff(iSession, 1) = shifts1.diff;
    
    
    
    
    
    
    
    

end % end of function