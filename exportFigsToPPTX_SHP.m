function exportFigsToPPTX_SHP(filename)
% Modified from Ingie Hong's ExportFigsToPPTX script
%
% Exports all open MATLAB figures to PPTX in vector graphic
% format.
% 
%
% Based on Example 1 of exportToPPTX
%
% Ingie Hong, Johns Hopkins Medical Institute, 2016

if nargin<1
    filename=datestr(now,'yymmdd_HHMMSS');
end

sortfigs % Sort figures in order generated

%% SAVE VECTOR FIGURES
%% Start new presentation
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end

exportToPPTX('new','Dimensions',[10 7.5], ...
    'Title',['MATLAB log on ' datestr(now,'yymmdd')], ...
    'Author','Soo Hyun Park', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','Automatically generated PPT file to save results figures');

% Additionally background color for all slides can be set as follows:
% exportToPPTX('new','BackgroundColor',[0.5 0.5 0.5]);


%% Add some slides
%figH = figure('Renderer','zbuffer'); mesh(peaks); view(0,0);
hfigs=get(0, 'Children');

for islide=1:length(hfigs)
    %figure(hfigs(islide))  
    slideNum = exportToPPTX('addslide');
    
    warning('off','MATLAB:print:CustomResizeFcnInPrint')
    saveas(hfigs(islide),'vectorFile','png'); %'emf');

    exportToPPTX('addpicture','vectorFile.png'); %'vectorFile.emf');
%     exportToPPTX('addtext',sprintf('Figure %d',islide));
end   
 
%% Save presentation and close presentation -- overwrite file if it already exists
% Filename automatically checked for proper extension
newFile = exportToPPTX('saveandclose', filename);
fprintf('New file has been saved: <a href="matlab:winopen(''%s'')">%s</a>\n',newFile,newFile);

%% Save MATLAB diary 
%diary([filedate '.log']) % has to be pre-initiated

delete('vectorFile.png');

clear islide
clear slideNum
clear isOpen
clear fileStats
clear newFile
clear hfigs

% Alternatively you can:
% exportToPPTX('save','example');
% exportToPPTX('close');


