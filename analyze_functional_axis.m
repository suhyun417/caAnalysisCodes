function result = analyze_functional_axis(x, y, PC, scale, nbins)
% Purpose: Evaluate periodic changes of PC scores along x or y axis of FOV%
% Input
%      x, y         - coordinates (in pixel)
%      PC            - scores of PCs (ex: PC1, PC2)
%      scale        - µm per pixel
%      nbins        - number of bins (default = 30)

if nargin < 5
    nbins = 30;
end

% convert pixels to micrometers
x_um = x * scale;
y_um = y * scale;

% result struct
result = struct();

% -------- axis loop --------
for axis = ["x", "y"]
    if axis == "x"
        coord = x_um;
    else
        coord = y_um;
    end

    % 1D binning
    edges = linspace(min(coord), max(coord), nbins+1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    bin_idx = discretize(coord, edges);
    binned_PC = accumarray(bin_idx, PC, [nbins 1], @mean, NaN);

    % get rid of NaNs
    valid = ~isnan(binned_PC);
    binned_PC = binned_PC(valid);
    centers = centers(valid);
    n = length(binned_PC);

    % FFT
    dx = mean(diff(centers));
    Y = fft(binned_PC - mean(binned_PC));
    freqs = (0:n-1) / (n*dx);
    power = abs(Y).^2;

    % dominant frequency (remove DC)
    [pmax, idx] = max(power(2:floor(n/2)));
    dom_freq = freqs(idx + 1);
    dom_lambda = 1 / dom_freq;

    % save
    result.(axis).freq = dom_freq;
    result.(axis).lambda = dom_lambda;
    result.(axis).power = pmax;
    result.(axis).binned_PC = binned_PC;
    result.(axis).centers = centers;
    result.(axis).power_spectrum = power(1:floor(n/2));
    result.(axis).freqs = freqs(1:floor(n/2));
end
end