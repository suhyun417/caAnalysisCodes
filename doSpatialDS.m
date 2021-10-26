function [Yf_resize, paramSpatialDS] = doSpatialDS(nameSubj, dateSession, fname, spatialFactor, method_DS)

% Perform spatial downsampling on .tif file 

addpath('/projects/parksh/_toolbox/TIFFstack'); % for loadtiff function

% set some default parameters
if nargin < 4
    spatialFactor = 4;
    method_DS = 'box'; %'bilinear'; %'bicubic'; %'nearest';
elseif nargin < 5
    method_DS = 'box';     
end

if str2num(dateSession) < 20191121
    dirRawdata = '/archive_rawdata1/parksh/calciumImaging/';
else
    dirRawdata = '/rawdata/parksh/calciumImaging/';
end
dirRawData_session =  fullfile(dirRawdata, [dateSession, '_', nameSubj]);

% load the file
fprintf(1, '   Perform spatial downsampling on %s...\n', fname);
Yf = loadtiff(fullfile(dirRawData_session, fname));
Yf = single(Yf);
[d1,d2,T] = size(Yf);

if rem(d1, spatialFactor) > 0 || rem(d2, spatialFactor) > 0
    Yf = Yf(1:floor(d1/spatialFactor).*spatialFactor, 1:floor(d2/spatialFactor).*spatialFactor, :);
    [d1,d2,T] = size(Yf);
end

% Spatial downsampling using imresize, frame-by-frame
Yf_resize = NaN(d1/spatialFactor, d2/spatialFactor, T);
for iFrame = 1:T
    Yf_resize(:,:,iFrame) = imresize(Yf(:,:,iFrame), 1/spatialFactor, method_DS);
end
fprintf(1, '     ...DONE! \n');


paramSpatialDS.spatialFactor = spatialFactor;
paramSpatialDS.imresize_method = method_DS;
paramSpatialDS.dim_DS = size(Yf_resize);
paramSpatialDS.dim_org = [d1, d2, T];

%        
%         clear Yf Yf_resize
%         
%         name_spatialDS = [nameRun, '_sDS_cat']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
%         save(fullfile(dirPreProcData_session, name_spatialDS), 'Yf_cat', 'paramPreproc');
%         fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled and concatenated image was saved as .mat\n', iRun, nRun, nameRun);
%         fastTiffStackWrite(fullfile(dirPreProcData_session, [name_spatialDS, '.tif']), Yf_cat);
% %         saveastiff(Yf_cat, fullfile(dirPreProcData_session, [name_spatialDS, '.tif']))
%         fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled and concatenated image was saved as .tif\n', iRun, nRun, nameRun);
%         
%     end