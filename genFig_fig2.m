% genFig_fig2.m
% 
% 2024/04/22 SHP
% Generate Figure 2 of marmoset calcium imaging manuscript
% a) example FOV and calcium traces
% b) longitudinal registration of the FOVs in two animals: extracted
% sources from each day, superimposed source boundaries over days
% c) entire population of cells (colored to indicate longitudinal
% tracking?)
% d) Histogram of tracking durations.  Need to think about whether total days, one for consecutive, or longest span.

directory = setDir_shp;
dirProjects = directory.dirProjects;
dirProcdata = directory.dirProcdata;
dirRawdata = directory.dirRawdata;