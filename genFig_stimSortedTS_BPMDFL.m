% genFig_stimSortedTS_BPMDFL.m
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


nameSubj = 'Tabla'; % 'Max'; %'Tabla';
% setDateSession = {'20191209'};

% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

clear infoSession


flagSavePPTX = 1; % 1 if you want to save figures to PPTX

for iSession = 2 %1:length(setDateSession)
    
    close all;
    
    dateSession = setDateSession{iSession}; %'20191205';
    
    dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
    switch lower(nameSubj)
        case 'tabla'
            dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
        case 'max'
            dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
    end
    
    % load the time series
    load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/BPM_ts_tML.mat', nameSubj, dateSession), 'tS_session_stim', 'stimTiming_BPM')
    load(sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Session/%s/DFL_ts_tML.mat', nameSubj, dateSession), 'tS_session', 'stimTiming_DFL')
    tS_session_DFL = tS_session;
    clear tS_session

    
    % load the source extraction data
    addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
    cnmfe_setup;
    d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
    
    %%
    condName_BPM = {'human face', 'marmoset face', 	'marmoset body', 'scene', 'non familiar object', 'hands and catcher',...
        'phase scrambled', 'space scrambled', 'grating', 'random dot motion'};
    
    setCond = [1 2 5 6 10]; %cat(2, 11:15, 21:25, 51
    nImage = 5;
    
    catCondMat = cat(1, stimTiming_BPM.condMat);
    curSetCond = unique(catCondMat(:,1));
    if sum(ismember(setCond, curSetCond))~=5
        setCond = curSetCond;
    end
    
    % temp = repmat([1:nImage], length(setCond), 1);
    % tempCondMat = cat(2, repmat(setCond', nImage, 1), temp(:));
    % setStimID = sort(tempCondMat(:,1)*10 + tempCondMat(:,2));
    setCondName = condName_BPM(setCond); %{infoTrial.infoStim([1:6:25]).nameCondition};
    
%     neuron_b = neuron.batches{1}.neuron;
%      [Cn, PNR] = neuron.correlation_pnr_parallel();
    thr = 0.3; % the lower the smaller (more centralized) the contour
    Coor = neuron.get_contours(thr);
    imgFOV = neuron.Cn.*neuron.PNR;
    
%     %         figure;
    neuron.show_contours([], [], imgFOV, 'true');
    
    figure;
    set(gcf, 'Color', 'w', 'Position', [583.5000  202.5000  753.0000  675.0000]) % Tabla [695   210   530   660]); %Max [583.5000  202.5000  753.0000  675.0000]) % Tabla
    [center] = neuron.estCenter();
    imagesc(imgFOV); colormap(gray);
    hold on
    plot(center(:,2), center(:, 1), 'r.');
    text(center(:,2)+1, center(:,1), num2str([1:length(neuron.ids)]'), 'Color', 'w');
    axis off
    
    
    % summary responses
    matAmpCellStim = reshape(cat(1, tS_session_stim.avgAmp), size(tS_session_stim));
    figure;
    set(gcf, 'Color', 'w', 'Position', [680         602        1080         376])
    imagesc(matAmpCellStim')
    colormap(jet)
    set(gca, 'CLim', [-1 1].*2)
    xlabel('Cells')
    ylabel('Stimulus')
    %         condName = {infoTrial.infoStim([1:6:25]).nameCondition};
    set(gca, 'YTickLabel', setCondName)
    colorbar;
    
     if flagSavePPTX % make a new ppt and save open figures so far
        % save figures
        addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
        addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
        
        fname_pptx = sprintf('%s_%s_visResp', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
        exportFigsToPPTX_SHP(fname_pptx);
     end
    
    cMap_run = bone(length(stimTiming_BPM)+1); % for each run
    cMap_cond = cool(length(setCond));
    catStimIDSession = cat(1, tS_session_stim(1,:).idStim);
    
    fig_BPMDFL = figure;
    set(fig_BPMDFL, 'Color', 'w', 'Position', [67          58        1144         917]);

    for iCell = 1:size(tS_session_stim,1)
        figure(fig_BPMDFL); clf;
        clear sp_stim
        
        for iCond = 1:length(setCond)
            idCond = setCond(iCond);
            for iM = 1:nImage
                indStim = find(catStimIDSession == 10*idCond + iM);
                
                figure(fig_BPMDFL);
                sp_stim(iCond, iM) = subplot(length(setCond)+1, nImage+1, (iCond-1)*(nImage+1)+iM);
                set(gca, 'ColorOrder', cMap_run(tS_session_stim(iCell, indStim).indTrial_org(:,1), :));
                hold on;
                plot([-1:0.1:3.5], tS_session_stim(iCell, indStim).matTS_norm, 'LineWidth', 1);
                set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
                
                avgTS(:,iM) = mean(tS_session_stim(iCell, indStim).matTS_norm, 2);
                
                if iM==1
                    ylabel(sp_stim(iCond, iM), setCondName{iCond}, 'Color', cMap_cond(iCond, :));
                end
            end
            
            figure(fig_BPMDFL);
            sp_stim(iCond, nImage+1) = subplot(length(setCond)+1, nImage+1, (iCond-1)*(nImage+1)+nImage+1);
            cla;
            set(gca, 'ColorOrder', jet(nImage).*0.8);
            hold on;
            plot([-1:0.1:3.5], avgTS)
            set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
            
            avgTS_cond(:, iCond) = mean(avgTS, 2);
        end
        
        sp_stim(length(setCond)+1, 1) = subplot(length(setCond)+1, nImage+1, length(setCond)*(nImage+1)+1);
        imagesc(imgFOV)
        colormap(sp_stim(length(setCond)+1, 1), gray)
        hold on
        plot(Coor{iCell}(1,:), Coor{iCell}(2,:), 'c.')
        axis off
        
        sp_stim(length(setCond)+1, nImage+1) = subplot(length(setCond)+1, nImage+1, (length(setCond)+1)*(nImage+1));
        set(gca, 'ColorOrder', cMap_cond);
        hold on;
        plot([-1:0.1:3.5], avgTS_cond)
        set(gca, 'XLim', [-1 3.5], 'XTick', -1:1:3)
        ylabel('norm. resp (z)')
        xlabel('Time from stim on (s)')
        
        title(sp_stim(1, 3), sprintf('Cell #%d/%d', iCell, size(tS_session_stim,1)))
        
        xlabel(sp_stim(5, 3), 'Time from stim on (s)')
        
        % Movie response
        clear sp_dfl
        figure(fig_BPMDFL);   
        sp_dfl(1) = subplot(length(setCond)+1, nImage+1, length(setCond)*(nImage+1)+[2 3]);
        plot(squeeze(tS_session_DFL(1).matTS_norm(:,iCell,:)))
        hold on;
        line(repmat(200:200:1000, 2, 1), repmat(get(gca, 'YLim'), 5, 1)', 'Color', ones(1, 3).*0.7);
        text(70:200:1100, ones(1, 6).*max(get(gca, 'YLim')), {'Sc' 'HF' 'Op' 'MB' 'Obj' 'MF'})
        axis tight
        title(sp_dfl(1), sprintf('Movie %s', tS_session_DFL(1).idStim))
        
        sp_dfl(2) = subplot(length(setCond)+1, nImage+1, length(setCond)*(nImage+1)+[4 5]);
        plot(squeeze(tS_session_DFL(2).matTS_norm(:,iCell,:)))
        hold on;
        line(repmat(200:200:1000, 2, 1), repmat(get(gca, 'YLim'), 5, 1)', 'Color', ones(1, 3).*0.7);
        text(70:200:1100, ones(1, 6).*max(get(gca, 'YLim')), {'Obj' 'Op' 'MF' 'Sc' 'MB' 'HF'})
        axis tight
        title(sp_dfl(2), sprintf('Movie %s', tS_session_DFL(2).idStim))
        
        set(sp_dfl, 'XTick', 0:200:1200, 'XTickLabel', 0:20:120)
        h = findobj(sp_dfl, 'type', 'axes');
        set([h(1:2).YLabel], 'string', 'Norm. resp (z)')
        set([h(1:2).XLabel], 'string', 'Time (s)')
        
        if flagSavePPTX
            % save figures
            addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
            addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));

            warning('off','MATLAB:print:CustomResizeFcnInPrint')
            saveas(fig_BPMDFL,'vectorFile','png'); 
            
            exportToPPTX('open', fname_pptx);
            slideNum = exportToPPTX('addslide');
            exportToPPTX('addpicture','vectorFile.png'); %'vectorFile.emf');
            exportToPPTX('saveandclose', fname_pptx);
            
            delete('vectorFile.png');
        end
        
    end
    
    if flagSavePPTX

%         fprintf('New file has been saved: <a href="matlab:winopen(''%s'')">%s</a>\n',newFile,newFile);
        
        switch lower(nameSubj)
            case 'tabla'
                dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
            case 'max'
                dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
        end
        movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
        
    end
%     
%     if flagSavePPTX
%         % save figures
%         addpath(fullfile(dirProjects, '/_toolbox/exportToPPTX/'));
%         addpath(fullfile(dirProjects, '/_toolbox/imagetools/'));
%         
%         fname_pptx = sprintf('%s_%s', nameSubj, dateSession); % fullfile('procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, 'FOV1', sprintf('%s.pptx', dateSession));
%         exportFigsToPPTX_SHP(fname_pptx);
%         
%         switch lower(nameSubj)
%             case 'tabla'
%                 dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Tabla/FOV1';
%             case 'max'
%                 dest = '/procdata/parksh/_marmoset/invivoCalciumImaging/Max/FOV3';
%         end
%         movefile(sprintf('./%s_%s*.pptx', nameSubj, dateSession), dest);
%     end
    
end
