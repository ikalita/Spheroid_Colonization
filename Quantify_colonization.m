%% Quantification of spheroid colonization using dual fluorescence labeling
%
% This script quantifies GFP and mCherry fluorescence intensity in 
% spheroids from microscopy Z-stack images.
%
% Workflow:
% 1) Quantify mean GFP and mCherry fluorescence intensity within the
%    spheroid mask for each Z-slice.
% 2) Subtract the corresponding background intensity measured from
%    non-infected control spheroids.
% 3) Calculate total integrated GFP and mCherry fluorescence intensity
%    for the spheroid.
% 4) Calculate the GFP/mCherry fluorescence ratio.
% 5) Save the results.
%
% The spheroid segmentation mask should be generated using the
% Segment_spheroid.m script.

clear;
close all;
clc;

%% Adjustable parameters

% Microscopy conditions
pixelSize = 1.38;    % Pixel size in microns
zStep = 10;          % Z-step size in microns

% Input image and segmentation mask
fileName = 'example.tif';
maskFile = 'example_seg.mat';

% Number of channels in the input TIFF
numChannels = 3;

% Channel assignments
GFPChannel = 1;
mCherryChannel = 2;
nuclearChannel = 3;

% Averaged non-infected control
backgroundFile = 'average_noninfected_control.mat';

%% Load image stack

info = imfinfo(fileName);
numPages = numel(info);

% The number of Z-slices to analyze
sizeZ = 16;

height = info(1).Height;
width = info(1).Width;

img = zeros(height, width, sizeZ, numChannels, 'uint16');

for z = 1:sizeZ
    for c = 1:numChannels
        
        page = (z - 1) * numChannels + c;

        img(:, :, z, c) = imread(fileName, page);
    end
end

%% Load segmentation mask

S = load(maskFile);

masksSample = S.masksSample;

%% Load averaged non-infected control

BG = load(backgroundFile);

%% Visual inspection of segmentation mask

zMask = round(sizeZ / 2);

figMask = figure('Name', 'Loaded segmentation mask');
figMask.WindowScrollWheelFcn = @scrollMaskFcn;

showMaskSlice();

%% Quantify mean fluorescence intensity per Z-slice

GFP_sample = nan(sizeZ, 1);
mCherry_sample = nan(sizeZ, 1);

GFP_corrected = nan(sizeZ, 1);
mCherry_corrected = nan(sizeZ, 1);

maskArea = nan(sizeZ, 1);

for z = 1:sizeZ

    % Get spheroid mask for the current Z-slice
    maskS = masksSample(:, :, z);

    % Calculate mask area in pixels
    maskArea(z) = nnz(maskS);

    % Convert fluorescence images to double
    GFPimg = double(img(:, :, z, GFPChannel));
    mCherryimg = double(img(:, :, z, mCherryChannel));

    if any(maskS(:))

        % Mean fluorescence intensity within the spheroid mask
        GFP_sample(z) = mean(GFPimg(maskS));
        mCherry_sample(z) = mean(mCherryimg(maskS));

        % Subtract corresponding non-infected control background
        GFP_corrected(z) = ...
            GFP_sample(z) - BG.GFP_control_mean(z);

        mCherry_corrected(z) = ...
            mCherry_sample(z) - BG.mCherry_control_mean(z);
    end
end

%% Calculate total integrated fluorescence intensity

% Integrated intensity for each channel is calculated as:
% mean corrected intensity per pixel x mask area,
% summed across all Z-slices.

GFP_spheroid = sum(GFP_corrected .* maskArea, 'omitnan');
mCherry_spheroid = sum(mCherry_corrected .* maskArea, 'omitnan');

% Calculate GFP/mCherry ratio
if mCherry_spheroid > 0
    GFP_mCherry_ratio = GFP_spheroid / mCherry_spheroid;
else
    GFP_mCherry_ratio = NaN;
    warning('Corrected mCherry intensity is <= 0; GFP/mCherry ratio was not calculated.');
end

%% Display results

fprintf('Corrected GFP intensity    = %.2f\n', GFP_spheroid);
fprintf('Corrected mCherry intensity = %.2f\n', mCherry_spheroid);
fprintf('GFP/mCherry ratio           = %.2f\n', GFP_mCherry_ratio);

%% Save results

[~, sampleName, ~] = fileparts(fileName);

sampleQuantFile = [sampleName '_quant.mat'];

save(sampleQuantFile, ...
    'GFP_sample', ...
    'mCherry_sample', ...
    'GFP_corrected', ...
    'mCherry_corrected', ...
    'maskArea', ...
    'GFP_spheroid', ...
    'mCherry_spheroid', ...
    'GFP_mCherry_ratio', ...
    'pixelSize', ...
    'zStep', ...
    'fileName', ...
    'maskFile');

%% Save summary results

resultsFile = 'results.txt';

% Create header if the summary file does not exist
if ~exist(resultsFile, 'file')

    fid = fopen(resultsFile, 'w');

    fprintf(fid, ...
        'Sample\tGFP_integrated\tmCherry_integrated\tGFP_mCherry_ratio\n');

    fclose(fid);
end

% Read existing summary file
fid = fopen(resultsFile, 'r');
existingText = fread(fid, '*char')';
fclose(fid);

lines = splitlines(string(existingText));

% Check whether the sample is already present
writeLine = true;

for i = 1:numel(lines)

    fields = split(lines(i), '\t');

    if ~isempty(fields) && strcmp(fields(1), sampleName)

        fprintf('Sample "%s" is already present in the summary file.\n', ...
            sampleName);

        writeLine = false;
        break;
    end
end

% Append result if the sample is not already present
if writeLine

    fid = fopen(resultsFile, 'a');

    fprintf(fid, '%s\t%.4f\t%.4f\t%.4f\n', ...
        sampleName, ...
        GFP_spheroid, ...
        mCherry_spheroid, ...
        GFP_mCherry_ratio);

    fclose(fid);

    fprintf('Results appended to %s\n', resultsFile);
end

%% Display functions

function scrollMaskFcn(~, event)

    zMask = evalin('base', 'zMask');
    sizeZ = evalin('base', 'sizeZ');

    % Scroll through Z-slices
    zMask = zMask - event.VerticalScrollCount;
    zMask = max(1, min(sizeZ, zMask));

    assignin('base', 'zMask', zMask);

    showMaskSlice();
end

function showMaskSlice()

    zMask = evalin('base', 'zMask');
    img = evalin('base', 'img');
    masksSample = evalin('base', 'masksSample');
    sizeZ = evalin('base', 'sizeZ');
    nuclearChannel = evalin('base', 'nuclearChannel');

    figure(gcf);
    clf;

    % Display nuclear fluorescence image
    imshow(img(:, :, zMask, nuclearChannel), []);
    hold on;

    % Spheroid mask
    BW = masksSample(:, :, zMask);

    % Extract spheroid boundary
    B = bwboundaries(BW);

    if ~isempty(B)

        nPoints = cellfun(@length, B);
        [~, idx] = max(nPoints);

        B = B{idx};

        % Overlay segmentation boundary
        plot(B(:, 2), B(:, 1), ...
            'r', 'LineWidth', 2);
    end

    title(sprintf('Spheroid segmentation: Z-slice %d / %d', ...
        zMask, sizeZ));
end