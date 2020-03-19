% runInitPreproc_CNMFE_dox.m
%
% Initial preprocessing (spatial downsampling and motion correction) before cell segmentation using CNMFe
% for doxycycline dataset
% modified from "doInitPreproc_CNMFe.m"
% 2020/03/18 SHP

clear all;
% gcp; % for parallel processingls
addpath('/projects/parksh/_toolbox/TIFFstack');
addpath('/projects/parksh/_toolbox/NoRMCorre/');
addpath('/projects/parksh/_toolbox/Fast_Tiff_Write/');
addpath('/projects/parksh/_toolbox/imagetools/');

%% Session info & optional parameters
setSubj ={'Tabla', 'Max'};

% Optional params for each subject
paramCrop(1).flagCrop = 0; 
paramCrop(1).coordsCrop = [20 320 1 270]; % for Tabla_FOV1 [Xstart, Xend, Ystart, Yend] for cropping
paramCrop(2).flagCrop = 0; %1;
paramCrop(2).coordsCrop = [41 336 17 261]; % %for Max_FOV2 [Xstart, Xend, Ystart, Yend] for cropping
flagBPF = 1;
flagSaveFigure = 1;
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

for iSubj = 1:length(setSubj)
    nameSubj = setSubj{iSubj};
    
    % get session info
    [infoSession, opts] = readInfoSession_dox(nameSubj);

    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    for iSession = 1:nSession
        dateSession = setDateSession{iSession};

        dirRawdata = '/rawdata/parksh/calciumImaging/';
        dirRawData_session =  fullfile(dirRawdata, [dateSession, '_', nameSubj]);
        dirProcData_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirProcData_session_preproc = fullfile(dirProcData_session, '_preproc');
        mkdir(dirProcData_session_preproc);
        
        % List of runs to process: get the info from the xls file
        locRun = find(((indRun==iSession+1).*(infoSession.(6)>0))>0);
        listRun = infoSession.(2)(locRun);
        
        %     d_raw = dir(fullfile(dirRawData_session, 'recording_*.raw'));
        %     [listRun{1:length(d_raw)}] = deal(d_raw.name);
        
        nRun = length(listRun);
        for iRun = 1:nRun
            
            nameRun = listRun{iRun}; %char(regexp(listRun{iRun}, '\d{6}(?=.raw)', 'match'));
            
            %% Concatenate separated .tif files from the same recording and perform Spatial downsampling
            name_spatialDS = [nameRun, '_sDS_cat'];
            
            if ~exist(fullfile(dirProcData_session_preproc, [name_spatialDS, '.tif']), 'file')
                
                % get the name(s) of .tif files for this run to concatenate
                d_tif = dir(fullfile(dirRawData_session, ['recording_', dateSession, '_', nameRun, '*.tif']));
                [~, ind] = sortrows({d_tif.name}', 'descend'); % to lineup the files in a right order for concatenation
                d_tif = d_tif(ind);
                listFiles = {d_tif.name}';
                
                % load each tif file and perform preprocessing
                nFile = length(listFiles);
                Yf_cat = [];
                for iFile = 1:nFile
                    
                    fprintf(1, ':: Run #%d/%d (%s): File #%d/%d: Loading %s from /rawdata to preprocess...\n', iRun, nRun, nameRun, iFile, nFile, listFiles{iFile});
                    
                    spatialFactor = 2; %4;
                    method_DS = 'box';
                    [Yf_resize, paramSpatialDS] = doSpatialDS(nameSubj, dateSession, listFiles{iFile}, spatialFactor, method_DS);
                    
                    Yf_cat = cat(3, Yf_cat, Yf_resize);
                    Yf_cat = single(Yf_cat); % downsampling changes it to double, so make single again
                end
                
                clear Yf Yf_resize
                
                paramPreproc.filename = listFiles;
                paramPreproc.paramSpatialDS = paramSpatialDS;
                
                %             name_spatialDS = [nameRun, '_sDS_cat']; %[nameSubj, '_', nameRun, '_RigidMotCorr']; 
                save(fullfile(dirProcData_session_preproc, name_spatialDS), 'Yf_cat', 'paramPreproc');
                fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled and concatenated image was saved as .mat\n', iRun, nRun, nameRun);
                fastTiffStackWrite(fullfile(dirProcData_session_preproc, [name_spatialDS, '.tif']), Yf_cat);
                fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled and concatenated image was saved as .tif\n', iRun, nRun, nameRun);
                
                clear Y*
                
                save(fullfile(dirProcData_session_preproc, 'paramPreproc.mat'), 'param*')
            else
                fprintf(1, ':: Run #%d/%d (%s): Spatially downsampled file exists as %s\n',iRun, nRun, nameRun, [name_spatialDS, '.tif'])
                clear Y*
            end
            
            %% Crop the image
            name_spatialDS_crop = [nameRun, '_sDS_cat_c'];
            
            if paramCrop(iSubj).flagCrop && ~exist(fullfile(dirProcData_session_preproc, [name_spatialDS_crop, '.tif']), 'file')
                coordsCrop = paramCrop(iSubj).coordsCrop;
                
                Y = loadtiff(fullfile(dirProcData_session_preproc, [name_spatialDS, '.tif']));
                Y_c = Y(coordsCrop(3):coordsCrop(4), coordsCrop(1):coordsCrop(2), :);
                fastTiffStackWrite(fullfile(dirProcData_session_preproc, [name_spatialDS_crop, '.tif']), Y_c);
                fprintf(1, ':: Run #%d/%d (%s): Images are cropped and saved as .tif\n', iRun, nRun, nameRun);
                clear Y Y_c
                
                save(fullfile(dirProcData_session_preproc, 'paramPreproc.mat'), 'param*')
            else
                fprintf(1, ':: Run #%d/%d (%s): Cropping is not necessary or the file already exists\n',iRun, nRun, nameRun)
%                 clear Y*
            end

            
            %% Gaussian bandpass filtering
            %             name_spatialDS_crop_bpf = [nameRun, '*_bpf'];
            
            if flagBPF && isempty(dir(fullfile(dirProcData_session_preproc, [nameRun, '*_bpf.tif']))) % regardless of cropping
                
                if paramCrop(iSubj).flagCrop
                    if exist(fullfile(dirProcData_session_preproc, [nameRun, '_sDS_cat_c.tif']), 'file')
                        fname_in =  fullfile(dirProcData_session_preproc, [nameRun, '_sDS_cat_c.tif']);
                    else
                        fprintf(1, ':: Run #%d/%d (%s): Error: Bandpass Filtering: the cropping needs to be done first.\n', iRun, nRun, nameRun);
                        return;
                    end
                else
                    if exist(fullfile(dirProcData_session_preproc, [nameRun, '_sDS_cat.tif']), 'file')
                        fname_in =  fullfile(dirProcData_session_preproc, [nameRun, '_sDS_cat.tif']);
                    else
                        fprintf(1, ':: Run #%d/%d (%s): Error: Bandpass Filtering: spatially downsampled file doesn''t exist.\n', iRun, nRun, nameRun);
                    end
                end
                
                % parameters for Gaussian bandpass filter (following Inscopix convention)
                f0 = 0.005; % lowpass filter cut-off frequency
                f1 = 0.5; % highpass filter cut-off frequency
                [Y_BP, paramBPF] = doGaussBPF(fname_in, f0, f1);
                
                [p, nameIn, ext] = fileparts(fname_in);
                fastTiffStackWrite(fullfile(dirProcData_session_preproc, [nameIn, '_bpf.tif']), single(Y_BP));
                fprintf(1, ':: Run #%d/%d (%s): Bandpass Filtered images are saved as .tif\n', iRun, nRun, nameRun);
                clear Y_BP
                
                save(fullfile(dirProcData_session_preproc, 'paramPreproc.mat'), 'param*')
            else
                fprintf(1, ':: Run #%d/%d (%s): Bandpass Filtered image file exists. \n',iRun, nRun, nameRun)
%                 clear Y*
            end
            
            %% Motion correction
            %% Rigid motion correction          

            % create template if it doesn't exist yet, using the first run
            if ~exist(fullfile(dirProcData_session_preproc, 'mc_template.tif'), 'file') 
                fprintf(1, '     :: Motion Correction: template file doesn''t exist: creating one now...\n')
                
                d = dir(fullfile(dirProcData_session_preproc, [listRun{1}, '*_bpf.tif'])); % use bandpass filtered data
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
                YY = YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:);
                
                clear Yf
                
                % set options
                [p, nameIn, ext] = fileparts(fname);
                options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false); %, ...
                %sprintf('_mc_bound%d.tif', bound)]));
                
                % read initial batch and compute template
                init_batch = options_r.init_batch; %100;
                add_value = options_r.add_value; %0;
                us_fac = options_r.us_fac;
                max_shift = options_r.max_shift;
                
                interval = ceil(T/2-init_batch/2+1):floor(T/2+init_batch/2);
                sizY = size(YY);
                T = sizY(end);
                nd = length(sizY)-1; % determine whether imaging is 2d or 3d
                sizY = sizY(1:nd);
                
                Y_temp = YY(:,:,interval);
                data_type = class(Y_temp);
                Y_temp = single(Y_temp);
                
                fprintf('Registering the first %i frames just to obtain a good template....',init_batch);
                template = median(Y_temp,nd+1)+add_value;
                fftTemp = fftn(template);
                for t = 1:size(Y_temp,nd+1)
                    [~,Greg] = dftregistration_min_max(fftTemp,fftn(Y_temp(:,:,t)), us_fac, -max_shift, max_shift, options_r.phase_flag);
                    M_temp = real(ifftn(Greg));
                    template = template*(t-1)/t + M_temp/t;
                end
                template = template + add_value;
                fprintf('..done. \n');
                
                paramMC.paramHPF = paramHPF;
                paramMC.options_r = options_r;
                
                % save the template
                save(fullfile(dirProcData_session_preproc, 'mc_template.mat'), 'paramMC', 'template')
                fastTiffStackWrite(fullfile(dirProcData_session_preproc, 'mc_template.tif'), single(template));
                fprintf(1, '     ......Done! Motion Correction: template file is saved as %s...\n', fullfile(dirProcData_session_preproc, 'mc_template.tif'))
                
                save(fullfile(dirProcData_session_preproc, 'paramPreproc.mat'), 'param*')
                clear Y Yf YY Y_temp template
            end
            
            if isempty(dir(fullfile(dirProcData_session_preproc, [nameRun, '*_mc.tif'])))
                
                fprintf(1, '    :: Motion Correction: template file is being loaded from %s\n', fullfile(dirProcData_session_preproc, 'mc_template.tif'))
                template_in = loadtiff(fullfile(dirProcData_session_preproc, 'mc_template.tif'));
                
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
                tic; [M1,shifts1,template1] = normcorre_batch(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r, template_in); toc % register filtered data
                
                % apply shifts to full dataset and save it as tif
                options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound, 'output_type', 'tif', 'tiff_filename', fullfile(dirProcData_session_preproc, [nameIn '_mc.tif']));
                tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc
                
                paramMC.paramHPF = paramHPF;
                
                [cY,mY,vY] = motion_metrics(YY(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
                [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
                [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
                %             [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);
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
                
            else
                fprintf(1, ':: Run #%d/%d (%s): Motion corrected file exists. \n',iRun, nRun, nameRun);
            end 
          
        end % run
        
    end % session
    
end % subject



