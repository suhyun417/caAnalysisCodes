% doInitPreProc_beforeMotionCorrection.m
%
% 2019/09/12 SHP
% Perform initial preprocessing before the motion correction step,
% including:
%         1. Crop the image to exclude the edge/blurry part of the field of view  
%         2. Run rigid-body motion correction using NoRMCorre toolbox to get the timepoints to cut out
%         3. Using the output from the step #2, remove the parts with large motion (above certain criterion)
%                   -- NoRMCorre toolbox by Eftychios A. Pnevmatikakis
%                        (https://github.com/flatironinstitute/NoRMCorre.git)
%                        especially the demo_1p.m in the toolbox
%         4. Save the preprocessed data in procdata/parksh/_marmoset/invivoCalciumImaging/SUBJECT/DATE/_preproc

clear all;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/NoRMCorre/')

%% Session info & optional parameters
dateSession = '20190917'; %'20180612'; %'20180924'; %'20180628'; %25'; 
nameSubj = 'Tabla'; %'Hoppy';
dirRawData_session = fullfile('/rawdata/parksh/calciumImaging/', [dateSession, '_', nameSubj]); %20180529_Hoppy';
dirProcData_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession, '_preproc');
mkdir(dirProcData_session);
% dirSaveFile = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj);
% mkdir(dirSaveFile, dateSession);
% dirProcData_session = fullfile(dirSaveFile, dateSession);

% Optional params
flagNoRMC = 0; % non-rigid body correction is optional
flagSaveFigure = 1;
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

% List of runs to process
tempd = dir(dirRawData_session);
[listRun{1:16}] = deal(tempd(5:20).name);  % manual coding



%%
nRun = length(listRun);
for iRun = 1:nRun
    
    % Set the filename
    nameRun = listRun{iRun}; 
    d = dir(fullfile(dirRawData_session, nameRun));
    filename = fullfile(dirRawData_session, nameRun, d(3).name);   
    
    saveFilename = strcat(nameRun, '_preproc');
    
    %% read data and convert to double
    Yf = read_file(filename);
    Yf = single(Yf);
    [d1,d2,T] = size(Yf);
    
    %% perform some sort of deblurring/high pass filtering
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
        Y = imfilter(Yf,psf,'symmetric');
        bound = 0;
    end
    
    %% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
    max_shift = 20; % 8; %10; %15;
%     bound = 40;
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift', max_shift,'iter',1,'correct_bidir',false);
    
    %% register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    figure; set(gcf, 'Position', [100 100 865 235]);
    plot(shifts_r);
    title(sprintf('max shift: %d', max_shift))
    
    %% Take out the part where there was big motion and boundary
    bound = 30;
    crit_maxmotion = 8; % roughly 10um in Mightex system
    validTP = find(sum(abs(shifts_r)>crit_maxmotion, 2)<1); 
    
    Y_preproc = Yf(bound/2+1:end-bound/2,bound/2+1:end-bound/2,1:300); %Yf(bound/2+1:end-bound/2,bound/2+1:end-bound/2,validTP);
    
    
    % test: do the motion correction on the preprocessed data
    Yf = single(Y_preproc);
    [d1,d2,T] = size(Yf);
    
    %% perform some sort of deblurring/high pass filtering
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
        Y = imfilter(Yf,psf,'symmetric');
        bound = 0;
    end
    
    %% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
    max_shift = 20; % 8; %10; %15;
%     bound = 40;
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift', max_shift,'iter',1,'correct_bidir',false);
    
    %% register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    figure; set(gcf, 'Position', [100 100 865 235]);
    plot(shifts_r);
    title(sprintf('%s: after clipping out shifts larger than %d pixels', nameRun, crit_maxmotion))
    
    % exclude boundaries due to high pass filtering effects
    tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to the preprocessed dataset
    


%     %% compute metrics 
% [cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
% [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
% 
% [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
% [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);
% 
% %% plot rigid shifts and metrics
% shifts_r = squeeze(cat(3,shifts1(:).shifts));
% figure;
%     subplot(311); plot(shifts_r);
%         title('Rigid shifts','fontsize',14,'fontweight','bold');
%         legend('y-shifts','x-shifts');
%     subplot(312); plot(1:T,cY,1:T,cM1);
%         title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
%         legend('raw','rigid');
%     subplot(313); plot(1:T,cYf,1:T,cM1f);
%         title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
%         legend('raw','rigid');
%     
      if flagSaveFigure
        print(fig_rigid, fullfile(dirFig, [nameSubj, '_', nameRun, '_RigidMotCorr']), '-depsc')
    end
    
    %% Save motion corrected files
    name_RigidMC = [nameRun, '_RigidMotCorr']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
    save(fullfile(dirProcData_session, name_RigidMC), 'Mr', 'shifts_r');     
    fprintf(1, 'Run #%d/%d: motion corrected images were saved as .mat\n', iRun, nRun);
    saveastiff(Mr, fullfile(dirProcData_session, name_RigidMC))
    fprintf(1, 'Run #%d/%d: motion corrected images were saved as .tiff\n', iRun, nRun);
    
end
    