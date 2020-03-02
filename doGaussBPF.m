function [Y_BP, paramBPF] = doGaussBPF(fname, f0, f1)

% Perform 2-D Gaussian Bandpass Filtering on imaging data, frame-by-frame
% following Inscopix data software parameters
%   Usage: [Yf_BP, paramBPF] = doGaussBPF(fname, f0, f1)
%   Input
%       fname: filename with path
%       f0: low-pass cut-off frequency (default: 0.005)
%       f1: high-pass cut-off frequency (default: 0.5)
%   Output
%       Y_BP: bandpass-filtered image in 3D
%       paramBPF: parameters used for BPF


addpath('/projects/parksh/_toolbox/TIFFstack'); % for loadtiff function

if nargin < 2
    f0 = 0.005; 
    f1 = 0.5;
elseif nargin <3
    f1 = 0.5;
end

% sigma
s0 = sqrt(2*log(2))/(2*pi*f0); % lowpass
s1 = sqrt(2*log(2))/(2*pi*f1); % highpass

% load data
fprintf(1, '   Perform gaussian bandpass filtering on %s...\n', fname);
Y = loadtiff(fname);
[d1,d2,T] = size(Y);

% Spatial downsampling using imresize, frame-by-frame
Y_BP = NaN(d1, d2, T);
for iFrame = 1:T
    Y_BP(:,:,iFrame) = imgaussfilt(Y(:,:,iFrame), s1) - imgaussfilt(Y(:,:,iFrame), s0); % highpass filtered image - lowpass filtered image
end
fprintf(1, '     ...DONE! \n');

paramBPF.fq_cutoff = [f0, f1];
paramBPF.sigma = [s0, s1];
paramBPF.fun = 'imgaussfilt';

