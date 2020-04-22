function [infoT, opts] = readInfoSession(nameSubj)

% Read a spreadsheet from a .xls file containing calcium imaging session information
% into a table, using "readtable" function
%   Usage: [infoT, opts] = readInfoSession(nameSubj);
%   Input
%       - nameSubj: full name of subject (string)
%   Output
%       - infoT: table containing the full contents
%       - opts: options for the table
% 2020/01/22 SHP


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

% get session info from spreadsheet
switch lower(nameSubj)
    case 'tabla'
        fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Tabla_FOV1.xlsx');
%         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Tabla_FOV3.xlsx');
    case 'max'
        fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Max_FOV3.xlsx');
%         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Max_FOV2.xlsx');
end

opts = detectImportOptions(fname_infoSession);
opts = setvartype(opts, {'Date', 'ImagingFilename', 'ExpName', 'MLFilename', 'stimulus', 'notes'}, 'char');
opts = setvartype(opts, 'flagPreproc', 'double');
infoT = readtable(fname_infoSession, opts);