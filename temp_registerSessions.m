% temp_registerSession.m
%
% 2021/09/28 SHP
% Registration across daily sessions using NoRMCorre
% For a given FOV, register daily sessions to the first day, using first
% run as a template

clear all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/nifvault/NIFVAULT/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

addpath('/nifvault/NIFVAULT/projects/parksh/_toolbox/TIFFstack');
addpath('/nifvault/NIFVAULT/projects/parksh/_toolbox/NoRMCorre/');
addpath('/nifvault/NIFVAULT/projects/parksh/_toolbox/Fast_Tiff_Write/');
addpath('/nifvault/NIFVAULT/projects/parksh/_toolbox/imagetools/');
% gcp; % for parallel processingls

%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

dirFig = '/nifvault/NIFVAULT/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

iSubj = 1; %1:length(setSubj)
nameSubj = setSubj{iSubj};

% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% Load the reference image (first session)
dirRefImage = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
    nameSubj, setDateSession{1});
imgRef = loadtiff(fullfile(dirRefImage, 'mc_template.tif'));

for iSession = 2:nSession
    % Load session image to align to the reference image
    dateSession = setDateSession{iSession};
    dirSessionImage = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
    nameSubj, dateSession);
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
    paramMC.options_r = options_r;
    
    % register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre(imgSession(bound/2+1:end-bound/2,bound/2+1:end-bound/2),options_r, imgRef); toc % register filtered data
    %      [M_final,shifts,template,options,col_shift] = normcorre(Y,options,template);
    % apply shifts and save it as desired format described in the options
    tic; Mrg = apply_shifts(imgSession,shifts1,options_r,bound/2,bound/2); toc    

%     %% now apply non-rigid motion correction
%     % non-rigid motion correction is likely to produce very similar results
%     % since there is no raster scanning effect in wide field imaging    
%     options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
%         'grid_size',[30,30]*2,'mot_uf',4,'correct_bidir',false, ...
%         'overlap_pre',32,'overlap_post',32,'max_shift',20, ...
%         'output_type', 'mat'); % 'tif', 'tiff_filename', './registrationTest_nrgm.tif');
%     
%     tic; [M2,shifts2,template2] = normcorre(imgSession(bound/2+1:end-bound/2,bound/2+1:end-bound/2),options_nr,imgRef); toc % register filtered data
% % apply shifts and save it as desired format described in the options
%     tic; Mnrg = apply_shifts(imgSession,shifts2,options_nr,bound/2,bound/2); toc    

% Check the registration
figure;
subplot(1, 2, 1);
imshowpair(imgRef, imgSession);
title('Before Registration')
subplot(1, 2, 2);
imshowpair(imgRef, Mrg);
title('Rigid Motion Regiration')
% subplot(1, 3, 3);
% imshowpair(imgRef, Mnrg);
% title('Non-rigid motion Registration')

        %% is it possible to apply shifts directly to neuron.A?
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        
        % REF IMAGE
        d_sources2D_ref = dir(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
            nameSubj, setDateSession{1}));
        load(fullfile(d_sources2D_ref(1).folder, d_sources2D_ref(1).name));        
        imgFOV_ref = neuron.Cn.*neuron.PNR;
        [center_ref] = neuron.estCenter();
                
        % SESSION IMAGE
        dirProcdata_session = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/',...
            nameSubj, dateSession);
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        clear neuron
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        imgFOV_ses = neuron.Cn.*neuron.PNR;
        imgFOV_ses_reg = apply_shifts(imgFOV_ses, shifts1,options_r,bound/2,bound/2);
        
        figure
        subplot(1,2,1);
        imshowpair(imgFOV_ref, imgFOV_ses)
        title('Before Registration')
        subplot(1,2,2);
        imshowpair(imgFOV_ref, imgFOV_ses_reg)
        title('After Registration')

        [center_ses] = neuron.estCenter();
        shifts = squeeze(shifts1.shifts);
        center_ses_reg = center_ses + shifts';
        

        figure
        sp(1) = subplot(1,2,1);
        plot(center_ref(:,2), center_ref(:, 1), 'r.')
        hold on
        plot(center_ses(:,2), center_ses(:, 1), 'b.')
        sp(2) = subplot(1,2,2);
        plot(center_ref(:,2), center_ref(:, 1), 'r.')
        hold on
        plot(center_ses_reg(:,2), center_ses_reg(:, 1), 'g.')
        set(sp, 'YDir', 'reverse')
        [d1,d2] = size(neuron.Cn);
        set(sp, 'XLim', [0 d2], 'YLim', [0 d1])
        title(sp(1), 'Cell center: before registration')
        title(sp(2), 'Cell center: after registration')



    thr = 0.5; % the lower the smaller (more centralized) the contour
cellColor = [1 1 1];
widthContour = 1;
[d1,d2] = size(neuron.Cn);

figure;
imagesc(zeros(d1, d2)); % background
colormap(gray);
caxis([0 0.1]);
hold on;

CC = cell(size(neuron.A, 2),1);
CR = cell(size(neuron.A, 2),2);
for i = 1:size(neuron.A ,2)
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
title(sprintf('%s: %s', nameSubj, dateSession))

%     [cY,mY,vY] = motion_metrics(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
%     [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
%     [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
    
%     paramHPF.gSig = 7;
%     paramHPF.gSiz = 17;
%     paramHPF.bound = 40;
%     paramHPF.imfilter = 'imfilter(Yf,psf,''symmetric'')';
    
%     Yf = loadtiff(fname);
%     [d1,d2,T] = size(Yf);
%     
%     % perform some sort of deblurring/high pass filtering
%     gSig = paramHPF.gSig; % 7;
%     gSiz = paramHPF.gSiz; % 17;
%     psf = fspecial('gaussian', round(2*gSiz), gSig);
%     ind_nonzero = (psf(:)>=max(psf(:,1)));
%     psf = psf-mean(psf(ind_nonzero));
%     psf(~ind_nonzero) = 0;   % only use pixels within the center disk
%     %Y = imfilter(Yf,psf,'same');
%     %bound = 2*ceil(gSiz/2);
%     YY = imfilter(Yf,psf,'symmetric');
%     bound = paramHPF.bound; %40; %0;
%     YY = YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:);
%     
%     clear Yf
%     
%     % set options
%     [p, nameIn, ext] = fileparts(fname);
%     options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false); %, ...
%     %sprintf('_mc_bound%d.tif', bound)]));
%     
%     % read initial batch and compute template
%     init_batch = options_r.init_batch; %100;
%     add_value = options_r.add_value; %0;
%     us_fac = options_r.us_fac;
%     max_shift = options_r.max_shift;
%     
%     interval = ceil(T/2-init_batch/2+1):floor(T/2+init_batch/2);
%     sizY = size(YY);
%     T = sizY(end);
%     nd = length(sizY)-1; % determine whether imaging is 2d or 3d
%     sizY = sizY(1:nd);
%     
%     Y_temp = YY(:,:,interval);
%     data_type = class(Y_temp);
%     Y_temp = single(Y_temp);
%     
%     fprintf('Registering the first %i frames just to obtain a good template....',init_batch);
%     template = median(Y_temp,nd+1)+add_value;
%     fftTemp = fftn(template);
%     for t = 1:size(Y_temp,nd+1)
%         [~,Greg] = dftregistration_min_max(fftTemp,fftn(Y_temp(:,:,t)), us_fac, -max_shift, max_shift, options_r.phase_flag);
%         M_temp = real(ifftn(Greg));
%         template = template*(t-1)/t + M_temp/t;
%     end
%     template = template + add_value;
%     fprintf('..done. \n');
                

    %% now apply non-rigid motion correction
    % non-rigid motion correction is likely to produce very similar results
    % since there is no raster scanning effect in wide field imaging
    
    options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
        'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
        'overlap_pre',32,'overlap_post',32,'max_shift',20);
    
    tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
    tic; Mpr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile

    
    %%
    imgRef = loadtiff(fullfile(dirProcData_session_preproc, 'mc_template.tif'));
    
    d = dir(fullfile(dirProcData_session_preproc, [nameRun, '*_bpf.tif'])); % use bandpass filtered data
    fname = fullfile(dirProcData_session_preproc, d.name);
    
    paramHPF.gSig = 7;
    paramHPF.gSiz = 17;
    paramHPF.bound = 40;
    paramHPF.imfilter = 'imfilter(Yf,psf,''symmetric'')';
    
    Yf = loadtiff(fname);
    [d1,d2,T] = size(Yf);
    
    % perform some sort of deblurring/high pass filtering
    gSig = paramHPF.gSig; % 7;
    gSiz = paramHPF.gSiz; % 17;
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    %Y = imfilter(Yf,psf,'same');
    %bound = 2*ceil(gSiz/2);
    YY = imfilter(Yf,psf,'symmetric');
    bound = paramHPF.bound; %40; %0;
    
    % set options
    [p, nameIn, ext] = fileparts(fname);
    clear options_r
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',T,'max_shift',20,'iter',1,'correct_bidir',false); %, ...
    paramMC.options_r = options_r;
    
    % register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre_batch(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r, imgRef); toc % register filtered data
    
    % apply shifts to full dataset and save it as tif
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound, 'output_type', 'tif', 'tiff_filename', fullfile(dirProcData_session_preproc, [nameIn '_mc.tif']));
    tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc
    
    paramMC.paramHPF = paramHPF;
    
    [cY,mY,vY] = motion_metrics(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
    [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
    [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);

    %
    %% plot rigid shifts and metrics
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    fig_rigid = figure(100);
    subplot(211); plot(shifts_r);
    title('Rigid shifts','fontsize',14,'fontweight','bold');
    legend('y-shifts','x-shifts');
    subplot(212); plot(1:T,cY,1:T,cM1);
    title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
    legend('raw','rigid');
    %             subplot(313); plot(1:T,cYf,1:T,cM1f);
    %             title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
    %             legend('raw','rigid');
    
    print(fig_rigid, fullfile(dirProcData_session_preproc, [nameIn '_mc']), '-depsc')
    
    clear Y* M*
    save(fullfile(dirProcData_session_preproc, 'paramPreproc.mat'), 'param*')