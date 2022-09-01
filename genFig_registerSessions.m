% genFig_registerSessions.m
%
% 2022/08/30 SHP
% - visualization of longitudinal registration results
% - same neuron's responses across trials over days

clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/PROJECTS/parksh';
        dirProcdata = '/Volumes/NIFVAULT/PROCDATA/parksh';
        dirRawdata = '/Volumes/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

nameSubj = 'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = 3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

% % Load the reference image (first session)
% dirRefImage = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/_preproc',...
%     nameSubj, setDateSession{1}));
% imgRef = loadtiff(fullfile(dirRefImage, 'mc_template.tif'));


%% load saved files
% cell-center info pooled across days
fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_stackedCenter.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_stack, 'stackCellCenter')

% cell quality info 
fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_cellQC, 'infoCells')

% translational shift across days
fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));  
load(fname_shifts, 'shifts')
    

%% Stacked cells across days
isCell = ~isnan(stackCellCenter);
catCell = sum(isCell, 3);
figure;
imagesc(catCell)
% catCell(8:14, 102:107)
% stackCellCenter(8:14, 102:107, :)

cellID = [86 75 96 69 55 61 63 nan]; %[6 10 13 10 12 8 6 nan]; %[109 93 113 nan nan 76 74 nan]; %[91 81 104 75 nan 64 nan nan]; %[69 53 68 nan 35 nan nan nan]; %[21 12 21 18 nan 12 12 16]; %[48 33 45 34 nan 32 33 32]; %[34 20 32 nan nan 28 23 25]; %[118 nan 28 nan nan 17 17 nan]; %[83 68 92 61 48 55 58 50]; %from Max FOV3
nameSubj = 'Max'; %'Tabla';
FOV_ID = 3; %1;
% [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);

tempMatTS1 = []; tempMatTS2 = []; nTrial = [];
for iS = 1:length(cellID)

    if isnan(cellID(iS))
        continue;
    end

    dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
    % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files

    dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);

    load(fullfile(dirProcdata_session, 'DFL_ts_tML'));

    figure(200);
    subplot(2,1,1)
    title('Mov 1')
    plot(squeeze(tS_session(1).matTS_norm(:, cellID(iS), :)))
    hold on

    tempMatTS1 = cat(1, tempMatTS1, squeeze(tS_session(1).matTS_norm(:, cellID(iS), :))');        
    nTrial = cat(1, nTrial, size(tempMatTS1, 1));

    subplot(2,1,2)
    title('Mov 2')
    plot(squeeze(tS_session(2).matTS_norm(:, cellID(iS), :)))
    hold on

    tempMatTS2 = cat(1, tempMatTS2, squeeze(tS_session(2).matTS_norm(:, cellID(iS), :))');
end

tempValidS = ~isnan(cellID);

figure;
set(gcf, 'color', 'w')
subplot(2,1,1)
imagesc(tempMatTS1);
set(gca, 'XTickLabel', 20:20:120, 'YTick', nTrial, 'YTickLabel', setDateSession(tempValidS),'TickDir', 'out', 'Box', 'off')
title('Mov 1')
colormap(hot)
subplot(2,1,2)
imagesc(tempMatTS2);
set(gca, 'XTickLabel', 20:20:120, 'YTick', nTrial, 'YTickLabel', setDateSession(tempValidS),'TickDir', 'out', 'Box', 'off')
title('Mov 2')
xlabel('Time (s)')
colormap(hot)

% %% playing
% for i = 1:length(infoCells)
%     indlow10 = [];
%     indlow10 = find(infoCells(i).pnrs<10);
%     
%     [sortpnrs, indsort] = sort(infoCells(i).pnrs, 'descend');
%     indtop10 = indsort(1:10);
%     
%     figure(100)
%     plot(infoCells(i).cellCenter(indtop10,2), infoCells(i).cellCenter(indtop10, 1), '.', 'MarkerSize', 15);
%     set(gca, 'YDir', 'reverse')
%     hold on;
%     input('')
% end

% %% plot movie tseries of potential same cell
% cellID = [22 16 14 17 19 18 nan nan nan 5 9 24]; % from Tabla FOV1
% nameSubj = 'Tabla'; %'Max'; %'Tabla';
% FOV_ID = 1; %3; %1;
% [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
% 
% [c, ia, indRun] = unique(infoSession.(1), 'sorted');
% setDateSession = c(2:end); % 1st one is always empty
% nSession = length(setDateSession);
% tempMatTS1 = []; tempMatTS2 = [];
% for iS = 1:length(cellID)
%     
%     if isnan(cellID(iS))
%         continue;
%     end
%         
%     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
%     % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
%     
%     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%     
%     load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
%     
%     figure(200);
%     subplot(2,1,1)
%     title('Mov 1')
%     plot(squeeze(tS_session(1).matTS_norm(:, cellID(iS), :)))
%     hold on
%     
%     tempMatTS1 = cat(1, tempMatTS1, squeeze(tS_session(1).matTS_norm(:, cellID(iS), :))');
%     
%     subplot(2,1,2)
%     title('Mov 2')
%     plot(squeeze(tS_session(2).matTS_norm(:, cellID(iS), :)))
%     hold on
%     
%     tempMatTS2 = cat(1, tempMatTS2, squeeze(tS_session(2).matTS_norm(:, cellID(iS), :))');
% end
    

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
    paramMC.options_r = options_r;
    
    % register using the high pass filtered data and apply shifts to original data
    tic; [M1,shifts1,template1] = normcorre(imgSession(bound/2+1:end-bound/2,bound/2+1:end-bound/2),options_r, imgRef); toc % register filtered data
    %      [M_final,shifts,template,options,col_shift] = normcorre(Y,options,template);
    % apply shifts and save it as desired format described in the options
%     tic; Mrg = apply_shifts(imgSession,shifts1,options_r,bound/2,bound/2); toc    

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

        %% Apply shifts directly to neuron.A: yes you can just add the shifts
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        
%         % REF IMAGE
%         d_sources2D_ref = dir(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
%             nameSubj, setDateSession{1})));
%         load(fullfile(d_sources2D_ref(1).folder, d_sources2D_ref(1).name));  
% %         tName = ls(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
% %             nameSubj, setDateSession{1}))); % because "dir" doesn't work 
% %         load(strtrim(tName));
%         imgFOV_ref = neuron.Cn.*neuron.PNR;
%         [center_ref] = neuron.estCenter();
                
        % SESSION IMAGE
        clear neuron
        dirProcdata_session = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/',...
            nameSubj, dateSession));
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));        
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
%         tName_ses = ls(fullfile(dirProcdata_session, 'Sources2D_all*'));
%         load(strtrim(tName_ses))
        imgFOV_ses = neuron.Cn.*neuron.PNR;
        imgFOV_ses_reg = apply_shifts(imgFOV_ses, shifts1,options_r,bound/2,bound/2);
        
%         figure
%         subplot(1,2,1);
%         imshowpair(imgFOV_ref, imgFOV_ses)
%         title('Before Registration')
%         subplot(1,2,2);
%         imshowpair(imgFOV_ref, imgFOV_ses_reg)
%         title('After Registration')

        [center_ses] = neuron.estCenter();
        shifts = squeeze(shifts1.shifts);
        center_ses_reg = center_ses + shifts';
        

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

%% Possible code
% 1. Function performing longitudinal registration
% - Load reference image
% - For each following sessions
%   - load session image
%   - perform the rigid-body registration
%   - save the shifts (and registration parameters and other outcomes)
%       - for each animal and each FOV (e.g. Tabla_FOV1_shifts.mat)
% 2. Another script to extract cells across daily sessions using the shifts
% - Maybe it's better to pool in the pixel space. imagine a grid in n x n pixel
% resolution (e.g. 3 pixels? I can choose the criterion) and make a list of
% grid in columnar way. Then gather cell IDs from the REF and SES 
% 3. Some kind of cell sorting/exclusion need to be done to exclude
% non-valid sources (e.g. not circular, too low SNR/PNR). 
%   Example (from Sources2D.m):
%       nA = sqrt(sum(obj.A.^2));
%       nr = length(nA);
%       K = size(obj.C, 1);
%       % order neurons based on its circularity
%       tmp_circularity = zeros(K,1);
%       for m=1:K
%           [w, r] = nnmf(obj.reshape(obj.A(:, m),2), 1);
%           ky = sum(w>max(w)*0.3);
%           kx = sum(r>max(r)*0.3);
%           tmp_circularity(m) = abs((kx-ky+0.5)/((kx+ky)^2));
%       end
%       [~, srt] = sort(tmp_circularity, 'ascend');
% 0. The cell exclusion can be done first before pooling cells across days.
% 0. So from now on we will always load the valid cell index computed from
% the above code #3 and apply that to the time series etc. Dealing with
% indices will be tricky and requires extra attention.

% ## Other things
% %% trim spatial components
%         function [ind_small] = trimSpatial(obj, thr, sz)
%             % remove small nonzero pixels
%             if nargin<2;    thr = 0.01; end
%             if nargin<3;    sz = 5; end
%             
%             se = strel('square', sz);
%             ind_small = false(size(obj.A, 2), 1);
%             for m=1:size(obj.A,2)
%                 ai = obj.A(:,m);
%                 ai_open = imopen(obj.reshape(ai,2), se);
%                 
%                 temp = full(ai_open>max(ai)*thr);
%                 l = bwlabel(obj.reshape(temp,2), 4);   % remove disconnected components
%                 [~, ind_max] = max(ai_open(:));
%                 
%                 ai(l(:)~=l(ind_max)) = 0;
%                 if sum(ai(:)>0) < obj.options.min_pixel %the ROI is too small
%                     ind_small(m) = true;
%                 end
%                 obj.A(:, m) = ai(:);
%             end
%             %             ind_small = find(ind_small);
%             %             obj.delete(ind_small);
%         end
% %% keep spatial shapes compact
%         function compactSpatial(obj)
%             for m=1:size(obj.A, 2)
%                 ai = obj.reshape(obj.A(:, m), 2);
%                 ai = circular_constraints(ai);
%                 obj.A(:, m) = ai(:);
%             end
%         end


%% Contour visualization
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