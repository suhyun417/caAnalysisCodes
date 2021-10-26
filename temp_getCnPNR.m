% genFig_stimSortedTS_BPM.m
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

nameSubj = 'Tabla';


[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
% setDateSession = {'20191209'};

for iSession = 1:length(setDateSession)
    
    dateSession = setDateSession{iSession}; %'20191205';
    
    dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    dirFig = fullfile(dirProjects, '/0Marmoset/Ca/_labNote/_figs/');
        
    % load the source extraction data
    addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
    cnmfe_setup;
    d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
    
    fprintf(1, 'Session %d/%d: %s: Computing Cn and PNR .....\n', iSession, length(setDateSession), dateSession)
    tic; [neuron.Cn, neuron.PNR] = neuron.correlation_pnr_parallel(); 
    
    save(fullfile(d_sources2D(1).folder, d_sources2D(1).name),  'neuron', 'paramCNMFE', '-v7.3');
    T= toc;
    fprintf(1, 'Session %d/%d: %s: ............................Done! Elapsed time: %d sec\n', iSession, length(setDateSession), dateSession, T)
    
end
