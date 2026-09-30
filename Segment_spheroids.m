%% Segmentation of spheroids from Z-stack images
% This script segments spheroid boundaries using nuclear labelling fluorescence.
%
% Workflow:
% 1) Segment the spheroid in each Z-slice.
% 2) Visually inspect and manually correct segmentation masks if needed.
% 3) Save the final segmentation masks.
%
% Input:
%   - Multi-channel TIFF Z-stack
%
% Output:
%   - MATLAB (.mat) file containing segmentation masks and boundaries

clear;
close all;
clc;

%% Microscopy parameters

pixelSize = 1.38;  % Pixel size in microns
zStep = 10;        % Z-step size in microns

%% Input and output settings

fileName = 'example.tif';       % CHANGE FILE NAME
outputFolder = 'segmentation_masks';

% Number of channels in the input TIFF
numChannels = 3;

% Channel containing the nuclear label used for segmentation
segmentationChannel = 3;

%% Read image metadata

info = imfinfo(fileName);
numPages = numel(info);

% Number of Z-slices to analyze
sizeZ = 16;

height = info(1).Height;
width = info(1).Width;

img = zeros(height, width, sizeZ, numChannels, 'uint16');

%% Import Z-stack

for z = 1:sizeZ
    for c = 1:numChannels

        % TIFF pages are ordered as:
        % Z1-Ch1, Z1-Ch2, Z1-Ch3, Z2-Ch1, Z2-Ch2, Z2-Ch3, ...
        page = (z - 1) * numChannels + c;

        img(:, :, z, c) = imread(fileName, page);
    end
end

%% Automatic segmentation of spheroid boundaries
%
% The spheroid is segmented using the nuclear fluorescence channel.
% If the segmented area decreases by more than 20% compared with the
% previous slice, the previous mask is used as a starting point to build 
% a new mask.
% The candidate mask with the highest difference between mean intensity
% inside and outside the mask is selected.

numZ = size(img, 3);

boundariesSample = cell(numZ, 1);
masksSample = false(height, width, numZ);
areas = nan(numZ, 1);

for z = 1:numZ

    % Extract and normalize the segmentation channel
    I = mat2gray(img(:, :, z, segmentationChannel));

    % Smooth image to reduce small-scale intensity variations
    Iblur = imgaussfilt(I, 10);

    % Threshold the image
    BW = imbinarize(Iblur, graythresh(Iblur));

    % Close small gaps in the spheroid boundary
    BW = imclose(BW, strel('disk', 55));

    % Fill holes inside the segmented spheroid
    BW = imfill(BW, 'holes');

    % Remove small objects
    BW = bwareaopen(BW, 500);

    % Keep the largest connected component
    BW = bwareafilt(BW, 1);

    areaCurrent = nnz(BW);

    %% Check for potential segmentation failure

    if z > 1

        areaRatio = areaCurrent / areas(z - 1);

        % If the segmented area decreases by >20%, use tracking mode
        if areaRatio < 0.8

            fprintf('Slice %d: tracking mode\n', z);

            BWprev = masksSample(:, :, z - 1);

            % Generate candidate masks by slightly expanding or
            % shrinking the mask from the previous Z-slice
            candidates = {
                BWprev
                imdilate(BWprev, strel('disk', 1))
                imdilate(BWprev, strel('disk', 2))
                imdilate(BWprev, strel('disk', 3))
                imdilate(BWprev, strel('disk', 4))
                imerode(BWprev, strel('disk', 1))
                imerode(BWprev, strel('disk', 2))
                imerode(BWprev, strel('disk', 3))
                };

            scores = nan(length(candidates), 1);

            % Score each candidate based on the intensity difference
            % between the inside and outside of the mask
            for k = 1:length(candidates)

                BWcand = candidates{k};

                insideMean = mean(I(BWcand));
                outsideMean = mean(I(~BWcand));

                scores(k) = insideMean - outsideMean;
            end

            % Select the candidate with the highest intensity contrast
            [~, bestIdx] = max(scores);

            BW = candidates{bestIdx};
            areaCurrent = nnz(BW);
        end
    end

    %% Store segmentation results

    areas(z) = areaCurrent;
    masksSample(:, :, z) = BW;

    B = bwboundaries(BW);

    if ~isempty(B)
        nPoints = cellfun(@length, B);
        [~, idx] = max(nPoints);
        boundariesSample{z} = B{idx};
    end
end

%% Visual inspection and manual correction

zSample = round(numZ / 2);

figSample = figure('Name', 'Spheroid segmentation');
figSample.WindowScrollWheelFcn = @scrollSampleFcn;
figSample.KeyPressFcn = @keyPressFcn;

showSampleSlice();

%% Save segmentation results

[~, name, ~] = fileparts(fileName);

% Create output folder if it does not already exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Output file
maskFile = fullfile(outputFolder, [name '_seg.mat']);

% Save segmentation masks, boundaries, and microscopy parameters
save(maskFile, ...
    'masksSample', ...
    'boundariesSample', ...
    'pixelSize', ...
    'zStep');

fprintf('Segmentation results saved to:\n%s\n', maskFile);

%% Callback functions

function scrollSampleFcn(~, event)

    zSample = evalin('base', 'zSample');
    numZ = evalin('base', 'numZ');

    % Scroll through Z-slices
    zSample = zSample - event.VerticalScrollCount;
    zSample = max(1, min(numZ, zSample));

    assignin('base', 'zSample', zSample);

    showSampleSlice();
end

%% Display function

function showSampleSlice()

    zSample = evalin('base', 'zSample');
    img = evalin('base', 'img');
    boundariesSample = evalin('base', 'boundariesSample');
    numZ = evalin('base', 'numZ');
    figSample = evalin('base', 'figSample');

    figure(figSample);
    clf;

    % Display the segmentation channel
    imshow(img(:, :, zSample, 3), []);
    hold on;

    % Overlay the spheroid boundary
    if ~isempty(boundariesSample{zSample})
        B = boundariesSample{zSample};
        plot(B(:, 2), B(:, 1), 'r', 'LineWidth', 2);
    end

    title(sprintf('Spheroid segmentation: Z-slice %d / %d', ...
        zSample, numZ));
end

%% Keyboard control

function keyPressFcn(~, event)

    switch event.Key

        case 'a'
            addToMask();

        case 'r'
            removeFromMask();
    end
end

%% Add region to mask

function addToMask()

    zSample = evalin('base', 'zSample');
    masksSample = evalin('base', 'masksSample');

    BW = masksSample(:, :, zSample);

    title('Draw region to ADD. Double-click when finished.');

    h = drawpolygon('Color', 'g');
    wait(h);

    BWadd = createMask(h);

    % Add the selected region to the current mask
    BW = BW | BWadd;
    BW = imfill(BW, 'holes');

    masksSample(:, :, zSample) = BW;

    assignin('base', 'masksSample', masksSample);

    updateBoundary(zSample);
    showSampleSlice();
end

%% Remove region from mask

function removeFromMask()

    zSample = evalin('base', 'zSample');
    masksSample = evalin('base', 'masksSample');

    BW = masksSample(:, :, zSample);

    title('Draw region to REMOVE. Double-click when finished.');

    h = drawpolygon('Color', 'y');
    wait(h);

    BWremove = createMask(h);

    % Remove the selected region from the current mask
    BW(BWremove) = false;

    masksSample(:, :, zSample) = BW;

    assignin('base', 'masksSample', masksSample);

    updateBoundary(zSample);
    showSampleSlice();
end

%% Update boundary after manual correction

function updateBoundary(zSample)

    masksSample = evalin('base', 'masksSample');
    boundariesSample = evalin('base', 'boundariesSample');

    BW = masksSample(:, :, zSample);
    B = bwboundaries(BW);

    if ~isempty(B)
        nPoints = cellfun(@length, B);
        [~, idx] = max(nPoints);
        boundariesSample{zSample} = B{idx};
    else
        boundariesSample{zSample} = [];
    end

    assignin('base', 'boundariesSample', boundariesSample);
end