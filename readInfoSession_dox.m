function [infoT, opts] = readInfoSession_dox(nameSubj)

% Read a spreadsheet from a .xls file containing calcium imaging session information
% into a table, using "readtable" function
% Modified from "readInfoSession.m"
%   Usage: [infoT, opts] = readInfoSession_dox(nameSubj);
%   Input
%       - nameSubj: full name of subject (string)
%   Output
%       - infoT: table containing the full contents
%       - opts: options for the table
% 2020/03/06 SHP


ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/ntis';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/ntis';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end

% get session info from spreadsheet
switch lower(nameSubj)
    case 'tabla'
        fname_infoSession = fullfile(dirProjects, 'infoSession_Tabla_FOV3_Dox.xlsx');
    case 'max'
        fname_infoSession = fullfile(dirProjects, 'infoSession_Max_FOV3_Dox.xlsx');
end

opts = detectImportOptions(fname_infoSession);
opts = setvartype(opts, {'Date', 'ImagingFilename', 'ExpName', 'MLFilename', 'stimulus', 'notes'}, 'char');
opts = setvartype(opts, 'flagPreproc', 'double');
infoT = readtable(fname_infoSession, opts);