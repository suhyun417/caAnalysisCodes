% doMotionCorrection.m
%
% 2018/05/23 SHP
% Perform motion correction on the 1p calcium imaging data
% using NoRMCorre toolbox by Eftychios A. Pnevmatikakis
% (https://github.com/flatironinstitute/NoRMCorre.git)
% especially the demo_1p.m in the toolbox
%   -- loads file from /rawdata/parksh/calciumImaging/SESSIONNAME/RUNNAME
%   -- perform motion correction: rigid-body correction (default) or
%   non-rigid-body correction (optional)
%   -- saves motion corrected file in
%   /procdata/parksh/_marmoset/invivoCalciumImaging/SUBJECT/
%   -- saves figures of before/after motion correction in
%   /projects/parksh/0Marmoset/Ca/_labNote/_figs

clear all;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/NoRMCorre/')

%% Session info & optional parameters
dateSession = '20190926'; % '20190910'; %'20180612'; %'20180924'; %'20180628'; %25'; 
nameSubj = 'Tabla'; %'Hoppy';
dirRawData_session = ['/rawdata/parksh/calciumImaging/', dateSession, '_', nameSubj]; %20180529_Hoppy';
dirSaveFile = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj);
mkdir(dirSaveFile, dateSession);
dirProcData_session = fullfile(dirSaveFile, dateSession);

% Optional params
flagNoRMC = 0; % non-rigid body correction is optional
flagSaveFigure = 1;
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

tempd = dir(dirRawData_session);
clear listRun
[listRun{1:7}] = deal(tempd(9:15).name);  % manual coding
% listRun = {'20180604_20x_100msec_3p0_blockdesign_run1_1';  
% '20180604_20x_100msec_3p0_blockdesign_run2_1';
% '20180604_20x_100msec_3p0_blockdesign_run3_1'; 
% '20180604_20x_100msec_3p0_blockdesign_run4_1';
% '20180604_20x_100msec_3p0_blockdesign_run5_1'; 
% '20180604_20x_100msec_3p0_blockdesign_run6_1';
% '20180604_20x_100msec_3p0_restingstate_run1_1';
% '20180604_20x_100msec_3p0_movieblock_run1_1';
% '20180604_20x_100msec_3p0_restingstate_run2_1';
% '20180604_20x_100msec_3p0_movieblock_run2_1';
% '20180604_20x_100msec_3p0_movieblock_run3_1';
% '20180604_20x_100msec_3p0_movieblock_run4_1';
% '20180604_20x_100msec_3p0_restingstate_run3_2';
% };



%%
nRun = length(listRun);
for iRun = 1:nRun
    nameRun = listRun{iRun}; % '20180524_20x_100msec_8p0_restingstate_green_1';
    % nameSubj = 'Hoppy';
    
    d = dir(fullfile(dirRawData_session, nameRun));
    filename = fullfile(dirRawData_session, nameRun, d(3).name);   
    
    
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
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
    
    %% register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
    tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
    % apply shifts on the whole movie
    
    %% compute metrics
    [cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
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
    
%     save(fullfile(dirProcData_session, name_RigidMC), 'Mr', 'shifts_r'); 
    
    if flagSaveFigure
        print(fig_rigid, fullfile(dirFig, [nameSubj, '_', nameRun, '_RigidMotCorr']), '-depsc')
    end
    
    %% Save motion corrected files
    name_RigidMC = [nameRun, '_RigidMotCorr']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; %'20180521_Hoppy_10x_100msec_5p_2_RigidMotCorr.tif';
    save(fullfile(dirProcData_session, name_RigidMC), 'Mr', 'shifts_r');     
    fprintf(1, 'Run #%d/%d: motion corrected images were saved as .mat\n', iRun, nRun);
    saveastiff(Mr, fullfile(dirProcData_session, name_RigidMC))
    fprintf(1, 'Run #%d/%d: motion corrected images were saved as .tiff\n', iRun, nRun);
    
    if flagNoRMC % non-rigid body correction is optional
        
        %% now apply non-rigid motion correction
        % non-rigid motion correction is likely to produce very similar results
        % since there is no raster scanning effect in wide field imaging
        options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
            'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
            'overlap_pre',32,'overlap_post',32,'max_shift',20);
        tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
        tic; Mpr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile
        
        %% compute metrics
        [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
        [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);
        
        %% plot shifts
        shifts_r = squeeze(cat(3,shifts1(:).shifts));
        shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
        shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
        shifts_x = squeeze(shifts_nr(:,2,:))';
        shifts_y = squeeze(shifts_nr(:,1,:))';
        patch_id = 1:size(shifts_x,2);
        str = strtrim(cellstr(int2str(patch_id.')));
        str = cellfun(@(x) ['patch # ',x],str,'un',0);
        figNonRigid = figure;
        ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients for filtered data','fontsize',14,'fontweight','bold')
        set(gca,'Xtick',[],'XLim',[0,T-3])
        ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
        set(gca,'Xtick',[])
        ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
        xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax1,ax2,ax3],'x')
        
        if flagSaveFigure
            print(figNonRigid, fullfile(dirFig, [nameSubj, '_', nameRun, '_NonRigidMotCorr']), '-depsc')
        end
    end
    
 end

