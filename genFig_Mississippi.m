% genFig_Mississippi.m
%
% Create files & figures for Mississippi's recording in 2020/01/08
% Mississippi is Jim Pickel's animal who got intravenous AAV injection with
% hSyn-GCaMP6f

clear all;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/TIFFstack');
addpath('/projects/parksh/_toolbox/NoRMCorre/');
addpath('/projects/parksh/_toolbox/Fast_Tiff_Write/');
addpath('/projects/parksh/_toolbox/imagetools/');

dirProcdata_session = '/procdata/parksh/_marmoset/invivoCalciumImaging/Mississippi';

%% Session info & optional parameters
nameSubj = 'Mississippi';
dateSession = '20200108';

flagSaveFigure = 1;
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

listRun = {'200umFromLensSurface.tif', '250umFromLensSurface.tif', '300umFromLensSurface.tif'};

%
nRun = length(listRun);
for iRun = 1:nRun
    
    clear Yf Yf_resize
    Yf = loadtiff(fullfile(dirProcdata_session, listRun{iRun}));
    Yf = single(Yf);
    [d1,d2,T] = size(Yf);
    
    % Spatial downsampling using imresize
    spatialFactor = 4;
    method_DS = 'box'; %'bilinear'; %'box'; %'bicubic'; %'nearest'; % 'bicubic';
    if rem(d1, spatialFactor) > 0 || rem(d2, spatialFactor) > 0
        Yf = Yf(1:floor(d1/spatialFactor).*spatialFactor, 1:floor(d2/spatialFactor).*spatialFactor, :);
        [d1,d2,T] = size(Yf);
    end
    
    Yf_resize = NaN(d1/spatialFactor, d2/spatialFactor, T);
    %             fprintf(1, ':: Run #%d/%d (%s): File #%d/%d: Spatial downsampling with spatial factor %2.2f, method %s ...',...
    %                 iRun, nRun, nameRun, iFile, nFile, 1/spatialFactor, method_DS);
    for iFrame = 1:T
        Yf_resize(:,:,iFrame) = imresize(Yf(:,:,iFrame), 1/spatialFactor, method_DS);
    end
    
    [~, nameRun, ext] = fileparts(listRun{iRun});
    name_spatialDS = [nameRun, '_sDS']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
    fastTiffStackWrite(fullfile(dirProcdata_session, [name_spatialDS, '.tif']), single(Yf_resize));
    fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled and concatenated image was saved as .tif\n', iRun, nRun, nameRun);
    
    
    %% Rigid motion correction
    clear Yf
    Yf = Yf_resize;
    clear Yf_resize
    
    [d1,d2,T] = size(Yf);
    
    % perform some sort of deblurring/high pass filtering
    if (0)
        hLarge = fspecial('average', 40);
        hSmall = fspecial('average', 2);
        for t = 1:T
            Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
        end
        %Ypc = Yf - Y;
        bound = size(hLarge,1);
    else
        gSig = 7;
        gSiz = 17;
        psf = fspecial('gaussian', round(2*gSiz), gSig);
        ind_nonzero = (psf(:)>=max(psf(:,1)));
        psf = psf-mean(psf(ind_nonzero));
        psf(~ind_nonzero) = 0;   % only use pixels within the center disk
        %Y = imfilter(Yf,psf,'same');
        %bound = 2*ceil(gSiz/2);
        YY = imfilter(Yf,psf,'symmetric');
        bound = 0;
    end
    
    %% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
    
    paramPreproc.rigidMotionCorrection.options = options_r;
    
    %% register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre_batch(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
    tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
    % apply shifts on the whole movie
    
    %% compute metrics
    [cY,mY,vY] = motion_metrics(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
    [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
    [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
    [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);
    
    %% plot rigid shifts and metrics
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    fig_rigid = figure;
    subplot(311); plot(shifts_r);
    title('Rigid shifts','fontsize',14,'fontweight','bold');
    legend('y-shifts','x-shifts');
    subplot(312); plot(1:T,cY,1:T,cM1);
    title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
    legend('raw','rigid');
    subplot(313); plot(1:T,cYf,1:T,cM1f);
    title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
    legend('raw','rigid');
    
    clear M1 shifts1 template1
    clear cY mY vY cYf mYf vYf cM1 mM1 vM1 cM1f mM1f vM1f
    
    if flagSaveFigure
        print(fig_rigid, fullfile(dirProcdata_session, [nameSubj, '_', nameRun, '_sDS_RigidMC']), '-depsc') %fullfile(dirFig, [nameSubj, '_', nameRun, '_sDS_cat_RigidMC']), '-depsc')
    end
    
    name_RigidMC = [nameRun, '_sDS_RigidMC']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
    fastTiffStackWrite(fullfile(dirProcdata_session, [name_RigidMC, '.tif']), single(Mr));
    fprintf(1, ':: Run #%d/%d (%s): Motion corrected image was saved as .tif\n', iRun, nRun, nameRun);
    
end

d = dir(fullfile(dirProcdata_session, '*_sDS_RigidMC.tif'));
for iRun = 1:nRun
    [p, nameRun, ext] = fileparts(d(iRun).name);
    saveFileName = fullfile(dirProcdata_session, [nameRun '_dFF.mat']);
    saveFileName_tif = fullfile(dirProcdata_session, [nameRun '_dFF.tif']);
    
    Mr = loadtiff(fullfile(dirProcdata_session, [nameRun, '.tif']));
    %         load(fullfile(dirPreProcData_session, d_mat(iRun).name), 'Mr');
    
    Y_avg = mean(Mr, 3);
    Y = (Mr - repmat(Y_avg, 1, 1, size(Mr,3)))./repmat(Y_avg, 1, 1, size(Mr,3));
    Y_std = std(Y, [], 3);
    Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
    
    save(saveFileName, 'Y_avg', 'Y', 'Y_std', 'Ysiz');
    fprintf(1, ':: Run #%d/%d: dFF of run %s from session %s_%s is saved as .mat file\n', iRun, nRun, dateSession, nameSubj, d(iRun).name)
    fastTiffStackWrite(saveFileName_tif, single(Y));
    fastTiffStackWrite(fullfile(dirProcdata_session, [nameRun '_dFF_std.tif']), Y_std);
    fprintf(1, ':: Run #%d/%d: dFF of run %s from session %s_%s is saved as .tif file\n', iRun, nRun, dateSession, nameSubj, d(iRun).name)
    
    clear Mr Y_avg Y_dFF Y_std
end

