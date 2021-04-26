% demoEyeDataBPM_20210226_SHP.m
%
% an example script to show how to read eye data stored in MonkeyLogic bhv2
% files, collected from marmoset calcium imaging experiments
% -- line-by-line execution is expected, not execution of an entire script

%% Define filename. 
% If you don't put in the full directory, Matlab will assume that this file is in the current directory
% (current directory  = where you're executing this line)
filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have

%% Read the file
data = mlread(filename);

% check the size of the entire struct, which corresponds to the number of trials 
length(data)
size(data, 2)
numTrial = length(data);


%% Extract all the conditions
% using "cat" function ("help cat" to see the syntax)
trial_conditions = cat(1, data.Condition); 

% condition info: sometimes you need to manually hard-code things like this..
% semi-colon inbetween items make this as a column vector (i.e.
% concatenation along the first dimension): see "help cat"
condName = {'human face'; 'marmoset face';	'marmoset body'; 'scene'; 'non familiar object'; 'familiar object';...
    'phase scrambled'; 'space scrambled'; 'grating'; 'random dot motion'};

% assessing data in cell array
condName{1} % this gives you the content "inside" of that cell
condName(1) % this gives you the cell itself


%% Which trials were "marmoset face" trials?
condition_id = 2; % for marmoset face (refer to the hard-coded "condName" above)
trial_id_marmosetFace = find(trial_conditions == condition_id);

% you can put number directly as below
trial_id_marmosetFace = find(trial_conditions == 2);

% You can also find the corresponding "id" for conditions by finding out
% certain characters or pattern (=sequence) of characters 
% For those, you use "strfind" or "contains" function. Check out their help
% pages

% E.g. two ways to get the "id" for any "face" conditions
contains(condName, 'face') 
strfind(condName, 'face') % this will need one more step to convert the output to the ID


%% Homework
% Try selecting different sets of trials based on the condition you want. 
% Use methods described above to get the ID for the conditions you want to
% select (e.g. all the face conditions, all the marmoset conditions, all the object conditions...).
% You can use different logical expressions that we talked about (e.g. >, ~= ).
 





