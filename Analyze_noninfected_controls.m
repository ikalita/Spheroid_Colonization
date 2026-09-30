%% Quantification of background fluorescence in non-infected controls
%
% This script quantifies the mean GFP and mCherry fluorescence intensity
% within a spheroid segmentation mask for each Z-slice of a non-infected
% control sample.
%
% Workflow:
% 1) Load a multi-channel TIFF Z-stack.
% 2) Load the corresponding spheroid segmentation masks.
% 3) Calculate mean GFP and mCherry fluorescence within the mask for
%    each Z-slice.
% 4) Save the results.

clear;
close all;
clc;

%% Adjustable parameters

% Microscopy conditions
pixelSize = 1.38;    % Pixel size in microns
zStep = 10;          % Z-step size in microns

% Input files
fileName = 'noninfected_example.tif';       % CHANGE FILE NAME
maskFile = 'noninfected_example_seg.mat';   % CHANGE FILE NAME

% Number of channels in the input TIFF
numChannels = 3;

% Channel assignments
GFPChannel = 1;
mCherryChannel = 2;
nuclearChannel = 3;

%% Load image stack

info = imfinfo(fileName);
numPages = numel(info);

% the number of Z-slices to analyze
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

% Check mask dimensions
if ~isequal(size(masksSample), [height, width, sizeZ])
    error('Mask dimensions do not match the image stack.');
end

%% Visual inspection of loaded segmentation masks

zMask = round(sizeZ / 2);

figMask = figure('Name', 'Loaded segmentation mask');
figMask.WindowScrollWheelFcn = @scrollMaskFcn;

showMaskSlice();

%% Quantify mean GFP and mCherry fluorescence

GFP_control = nan(sizeZ, 1);
mCherry_control = nan(sizeZ, 1);
maskArea = nan(sizeZ, 1);

for z = 1:sizeZ

    % Get segmentation mask for the current Z-slice
    mask = masksSample(:, :, z);

    % Calculate mask area in pixels
    maskArea(z) = nnz(mask);

    % Convert image data to double for intensity calculations
    GFPimg = double(img(:, :, z, GFPChannel));
    mCherryimg = double(img(:, :, z, mCherryChannel));

    % Calculate mean fluorescence intensity within the mask
    if any(mask(:))
        GFP_control(z) = mean(GFPimg(mask));
        mCherry_control(z) = mean(mCherryimg(mask));
    end
end

%% Save quantification results

[~, name, ~] = fileparts(fileName);

outputFolder = 'non_infected_control/quant_results';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

quantFile = fullfile(outputFolder, [name '_quant.mat']);

save(quantFile, ...
    'GFP_control', ...
    'mCherry_control', ...
    'maskArea', ...
    'pixelSize', ...
    'zStep', ...
    'fileName');

fprintf('Quantification results saved to:\n%s\n', quantFile);

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

    % Segmentation mask
    BW = masksSample(:, :, zMask);

    % Extract spheroid boundary
    B = bwboundaries(BW);

    if ~isempty(B)

        % Select the largest boundary
        nPoints = cellfun(@length, B);
        [~, idx] = max(nPoints);

        B = B{idx};

        % Overlay segmentation boundary
        plot(B(:, 2), B(:, 1), ...
            'r', 'LineWidth', 2);
    end

    title(sprintf('Non-infected control: Z-slice %d / %d', ...
        zMask, sizeZ));
end