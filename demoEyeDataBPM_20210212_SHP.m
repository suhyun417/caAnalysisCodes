% demoEyeDataBPM_20210212_SHP.m
%
% an example code to show how to read eye data stored in MonkeyLogic bhv2
% files, collected from marmoset calcium imaging experiments

%% Define filename. 
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have

%% Read the file
data = mlread(filename);

%% Look at one trial
iTrial = 20; %which trial?
data(iTrial).AnalogData % check the fields within this struct

size(data(iTrial).AnalogData.Eye) % always a good idea to check the size of data you're assessing 
xData = data(iTrial).AnalogData.Eye(:,1); % retrieve first column and assign it to a variable named "xData"
yData = data(iTrial).AnalogData.Eye(:,2); % retrieve second column and assign it to a variable named "yData"
pData = data(iTrial).AnalogData.General.Gen1; % pupil data

%% Basic plot
figure;
plot(xData)

help plot

figure;
plot(xData, yData);

figure;
plot([xData, yData]);

figure;
plot(xData);
hold on;
plot(yData);


%% Homework
% H1. How many trials are there? How can you get that information from BPM data ?
% (hint: "size", or "length" function)
%
% H2. How can you extract all the conditions from all the trials?
% (hint: see how this line works: "aa_trialNumbers = cat(1, data.Trial);")
%
% H3. "data(iTrial).BehavioralCodes" contains code numbers and corresponding
% timing information. Our code number 20 means stimulus on, 55 means
% stimulus off. How can you get the timing of stimulus on and stimulus
% off from a given trial? Think about how you can extract the pupil data *during* the stimulus presentation
% (i.e. between stimulus on and off) (We will go over this next time)


