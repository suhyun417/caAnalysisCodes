% genFig_stimSortedTS_DFL.m
%

clear all; close all;
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

addpath('/projects/parksh/_toolbox')
addpath('/projects/parksh/_toolbox/imagetools')
addpath('/projects/parksh/_toolbox/exportToPPTX/')

flagSavePPTX = 1; % 1 if you want to save figures to PPTX


setNameSubj = {'Tabla', 'Max'};

for iSubj = 1:length(setNameSubj)
    
    nameSubj = setNameSubj{iSubj}; %'Tabla'; %'Max';
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    for iSession = 1:length(setDateSession)
        
        close all;
        
        dateSession = setDateSession{iSession}; %'20191205';
        
        dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
        
        % load the time series
        load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession), 'tS_session', 'stimTiming_DFL')
        
        % load the source extraction data
        addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        %%
        neuron_b = neuron.batches{1}.neuron;
        thr = 0.3; % the lower the smaller (more centralized) the contour
        Coor = neuron_b.get_contours(thr);
        imgFOV = neuron_b.Cn.*neuron_b.PNR;
        
        neuron_b.show_contours([], [], imgFOV, 'true');
        
        
        for iCell = 1:size(tS_session(1).matTS, 2) %length(setCell) %1:size(matTS_movie, 1)
            figCheck = figure;
            set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1380 400])
            
            sp(1) = subplot(2, 3, [1 4]);
            imagesc(imgFOV)
            colormap(sp(1), gray)
            hold on
            plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'c.')
            axis off
            
            sp(2) = subplot(2, 3, [2 3]);
            plot(squeeze(tS_session(1).matTS_norm(:,iCell,:)))
            hold on;
            %         plot(tS_session(1).avgTS_norm(:,iCell), 'k-')
            legend;
            axis tight
            title(sp(2), sprintf('Movie %s', tS_session(1).idStim))
            
            sp(3) = subplot(2, 3, [5 6]);
            plot(squeeze(tS_session(2).matTS_norm(:,iCell,:)))
            hold on;
            %         plot(tS_session(2).avgTS_norm(:,iCell), 'k-')
            legend;
            axis tight
            title(sp(3), sprintf('Movie %s', tS_session(2).idStim))
            
            set(sp(2:3), 'XTick', 0:200:1200, 'XTickLabel', 0:20:120)
            h = findobj(gcf, 'type', 'axes');
            set([h(1:2).YLabel], 'string', 'Norm. resp (z)')
            set([h(1:2).XLabel], 'string', 'Time (s)')
            
            p = mtit(sprintf('%s %s: Cell #%d/%d', dateSession, nameSubj, iCell, size(tS_session(1).avgTS, 2)));
            
            %         imagesc(zscore(squeeze(matTS_movie(iCell, :, :)))')
            %         colormap(hot)
            %         ylabel(SP(1), 'Viewing')
            %         set(SP(1), 'CLim', [0 8]);
            %
            %         F(iCell) = getframe(gcf);
            %         drawnow
        end
        
        % save figures
        if flagSavePPTX
            
            fname_pptx = sprintf('%s_%s_visResp', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
%             exportFigsToPPTX
            exportFigsToPPTX(fname_pptx);
            
            switch lower(nameSubj)
                case 'tabla'
                    dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
                case 'max'
                    dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
            end
            movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
        end
        
    end
end



%         fname_movie = sprintf('%s_%s_DFL_%dCells_tSeries_movie%s', dateSession, nameSubj, length(setCell), setMovie{iMovie});
%         writerObj = VideoWriter(fullfile(dirProcdata_session, fname_movie), 'MPEG-4');
%         writerObj.FrameRate = 1;
%         open(writerObj);
%         for i = 1:length(F)
%             frame = F(i);
%             writeVideo(writerObj, frame);
%         end
%         close(writerObj);


% 
% if flagSavePPTX
%     % save figures
%     addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
%     addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
%     
%     fname_pptx = sprintf('%s_%s', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
%     exportFigsToPPTX(fname_pptx);
%     
%     switch lower(nameSubj)
%         case 'tabla'
%             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
%         case 'max'
%             dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
%     end
%     movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
% end
% 
% end
