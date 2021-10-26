% analCa_Dox_neuron.m
%
% 2020/04/24 SHP

clear all; 

%% directory
dirProjects = '/projects/parksh';
dirProcdata = '/procdata/parksh';
dirRawdata = '/rawdata/parksh';

flagSaveFile = 1;

setSubj ={'Tabla', 'Max'};

%%
for iSubj = 1:length(setSubj)
    nameSubj = setSubj{iSubj};
    % get session info
    [infoSession, opts] = readInfoSession_dox(nameSubj);
    
    [c, ia, indRun] = unique(infoSession.(1), 'sorted');
    setDateSession = c(2:end); % 1st one is always empty
    nSession = length(setDateSession);
    
    clear infoSession
    
    % dox directory for each subject
    dirDox = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Dox', nameSubj);
    
    for iSession = 1:length(setDateSession)
        
        dateSession = setDateSession{iSession}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
        % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
        
        dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
        dirPreproc = fullfile(dirProcdata_session, '_preproc');
        
        dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');        
        
        if iSubj==2 && iSession == length(setDateSession) % for Max's last session
            resultsDox_neuron(iSession).dateSession = dateSession;
            resultsDox_neuron(iSession).nCell = 0;
            resultsDox_neuron(iSession).C_raw = 0;
            resultsDox_neuron(iSession).C = 0;
            resultsDox_neuron(iSession).A = 0;
            resultsDox_neuron(iSession).Cn = 0;
            resultsDox_neuron(iSession).PNR = 0;
            resultsDox_neuron(iSession).Coor = 0;
            
            resultsDox_neuron_paramCNMFE(iSession).critPNR = paramCNMFE.critPNR;
            resultsDox_neuron_paramCNMFE(iSession).sortPNR = paramCNMFE.sortPNR;
            resultsDox_neuron_paramCNMFE(iSession).min_pnr = paramCNMFE.min_pnr;
            resultsDox_neuron_paramCNMFE(iSession).critCn = paramCNMFE.critCn;
            resultsDox_neuron_paramCNMFE(iSession).sortCn = paramCNMFE.sortCn;
            resultsDox_neuron_paramCNMFE(iSession).min_corr = paramCNMFE.min_corr;
        else
        %% Load the source extraction data
        addpath(fullfile(dirProjects, '_toolbox/CNMF_E/'));
        cnmfe_setup;
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        if isempty(d_sources2D)
            copyfile(fullfile(dirProcdata_session, '_preproc', 'ConcatRuns_all_source_extraction/Sources2D_all*.mat'), dirProcdata_session)
        end
        d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
        load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
        
        resultsDox_neuron(iSession).dateSession = dateSession;
        resultsDox_neuron(iSession).nCell = size(neuron.A, 2);
        resultsDox_neuron(iSession).C_raw = neuron.C_raw;
        resultsDox_neuron(iSession).C = neuron.C;
        resultsDox_neuron(iSession).A = neuron.A;
        resultsDox_neuron(iSession).Cn = neuron.Cn;
        resultsDox_neuron(iSession).PNR = neuron.PNR;
        resultsDox_neuron(iSession).Coor = neuron.Coor;
        
        resultsDox_neuron_paramCNMFE(iSession).critPNR = paramCNMFE.critPNR;
        resultsDox_neuron_paramCNMFE(iSession).sortPNR = paramCNMFE.sortPNR;
        resultsDox_neuron_paramCNMFE(iSession).min_pnr = paramCNMFE.min_pnr;
        resultsDox_neuron_paramCNMFE(iSession).critCn = paramCNMFE.critCn;
        resultsDox_neuron_paramCNMFE(iSession).sortCn = paramCNMFE.sortCn;
        resultsDox_neuron_paramCNMFE(iSession).min_corr = paramCNMFE.min_corr;
        end
        
        if flagSaveFile
            fname_save = sprintf('doxResults_neuron_%s.mat', nameSubj);
%             save(fullfile(dirDox, fname_save), 'resultsDox_neuron');
            save(fullfile(dirDox, fname_save), 'resultsDox*');
        end
    
    end
    
end
        
        
        