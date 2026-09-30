%% Average non-infected control profiles
%
% This script averages GFP and mCherry fluorescence profiles from
% multiple non-infected control samples.

clear;
close all;
clc;

%% Microscopy conditions

pixelSize = 1.38;  % Pixel size in microns
zStep = 10;        % Z-step size in microns

%% Input folder

% Folder containing individual control quantification results
quantFolder = 'non_infected_control/quant_results';

% Find all quantification files
files = dir(fullfile(quantFolder, '*_quant.mat'));

nFiles = numel(files);

if nFiles == 0
    error('No quantification files found in %s.', quantFolder);
end

if nFiles ~= 3
    warning('Expected 3 control files, found %d.', nFiles);
end

%% Load first file to determine the number of Z-slices

firstFile = fullfile(files(1).folder, files(1).name);
S = load(firstFile);

numZ = numel(S.GFP_control);

GFP_all = nan(numZ, nFiles);
mCherry_all = nan(numZ, nFiles);

%% Load all control profiles

for i = 1:nFiles

    filePath = fullfile(files(i).folder, files(i).name);
    S = load(filePath);

    % Check that all controls have the same number of Z-slices
    if numel(S.GFP_control) ~= numZ || ...
            numel(S.mCherry_control) ~= numZ

        error('Different numbers of Z-slices detected in %s.', ...
            files(i).name);
    end

    GFP_all(:, i) = S.GFP_control;
    mCherry_all(:, i) = S.mCherry_control;
end

%% Calculate average control profiles

GFP_control_mean = mean(GFP_all, 2, 'omitnan');
mCherry_control_mean = mean(mCherry_all, 2, 'omitnan');

GFP_control_std = std(GFP_all, 0, 2, 'omitnan');
mCherry_control_std = std(mCherry_all, 0, 2, 'omitnan');

%% Save results

saveFile = fullfile(quantFolder, ...
    'average_noninfected_control.mat');

% Store names of the files included in the average
controlFiles = {files.name};

save(saveFile, ...
    'GFP_control_mean', ...
    'mCherry_control_mean', ...
    'GFP_control_std', ...
    'mCherry_control_std', ...
    'controlFiles', ...
    'pixelSize', ...
    'zStep');

fprintf('Average control profile saved to:\n%s\n', saveFile);