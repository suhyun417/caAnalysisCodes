function [] = doMotionCorrection_lowRAM(fname, flagSaveFigure)
%
% 2020/1/1 SHP
% Perform motion correction on the 1p calcium imaging data
% using NoRMCorre toolbox by Eftychios A. Pnevmatikakis
% (https://github.com/flatironinstitute/NoRMCorre.git)
% especially the demo_1p_low_RAM.m in the toolbox
%   -- saves motion corrected file and results figures in /_preproc under
%   /procdata/parksh/_marmoset/invivoCalciumImaging/SUBJECT/


% clear;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/NoRMCorre/')

%% read data and convert to double
frame = read_file(fname,1,1);
[d1,d2] = size(frame);

%% perform some sort of deblurring/high pass filtering
% The function does not load the whole file in memory. Instead it loads
% chunks of the file and then saves the high pass filtered version in a
% h5 file.

gSig = 7;
gSiz = 3*gSig;
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk

[filepath,file_name,ext] = fileparts(fname);
h5_name = fullfile(filepath,[file_name,'_filtered_data.h5']);
chunksize = 750;    % read 500 frames at a time
cnt = 1;
while (1)  % read filter and save file in chunks
    Yf = single(read_file(fname,cnt,chunksize));
    if isempty(Yf)
        break
    else
        YY = imfilter(Yf,psf,'symmetric');
        saveash5(YY,h5_name);
        cnt = cnt + size(YY,ndims(YY));
    end
    disp(cnt)
end
%% first try out rigid motion correction
% exclude boundaries due to high pass filtering effects
options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);

%% register using the high pass filtered data and apply shifts to original data
tic; [M1,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
% exclude boundaries due to high pass filtering effects

% if you save the file directly in memory make sure you save it with a
% name that does not exist. Change options_r.tiff_filename
% or options_r.h5_filename accordingly.

fname_out = fullfile(filepath, [file_name '_RigidMC']);
options_r = NoRMCorreSetParms('output_type', 'tif', 'tiff_filename', [fname_out '.tif']);

tic; Mr = apply_shifts(fname,shifts1,options_r); toc % apply shifts to full dataset and save as tif

% Y = single(Mr);
% Ysiz  = size(Y)';
shifts_r = squeeze(cat(3,shifts1(:).shifts));
% fprintf(1, ':: Saving motion corrected images as .mat file...\n')
% save([fname_out '.mat'], 'Y', 'Ysiz', 'shifts_r', 'options_r', '-v7.3')
fprintf(1, '.....Done!\n')

% you can only save the motion corrected file directly in memory by
% setting options_r.output_type = 'tiff' or 'h5' and selecting an
% appropriate name through options_r.tiff_filename or options_r.h5_filename

% % compute metrics
% [cY,mY,vY] = motion_metrics(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
% [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
% 
% [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
% [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

%% plot rigid shifts and metrics
shifts_r = squeeze(cat(3,shifts1(:).shifts));
fig_rigid = figure;
% subplot(311); 
plot(shifts_r);
title('Rigid shifts','fontsize',14,'fontweight','bold');
legend('y-shifts','x-shifts');
% subplot(312); plot(1:T,cY,1:T,cM1);
% title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
% legend('raw','rigid');
% subplot(313); plot(1:T,cYf,1:T,cM1f);
% title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
% legend('raw','rigid');

if flagSaveFigure
    print(fig_rigid, fname_out, '-depsc') %fullfile(dirFig, [nameSubj, '_', nameRun, '_sDS_cat_RigidMC']), '-depsc')
end

% name_RigidMC = [nameRun, '_sDS_cat_RigidMC']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
% save(fullfile(dirPreProcData_session, name_RigidMC), 'Mr', 'paramPreproc', 'shifts_r');
% fprintf(1, ':: Run #%d/%d (%s): Motion corrected image was saved as .mat\n', iRun, nRun, nameRun);
% saveastiff(Mr, fullfile(dirPreProcData_session, [name_RigidMC, '.tif']))
% fprintf(1, ':: Run #%d/%d (%s): Motion corrected image was saved as .tif\n', iRun, nRun, nameRun);


end

