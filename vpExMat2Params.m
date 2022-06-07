% Francois Kroll 17/04/2022

% small script to extract parameter data from .mat file created by
% Vp_Extract.m

% alternative would be to read the .mat file in R but this seems to be
% complicated

% Note, run it once per .mat file

%% select .mat file

% let user select file
[filename, pathname] = uigetfile('*.mat', 'Select .mat file','MultiSelect','off'); %Select files
if isequal(filename,0) %If no file is selected
    error('No Files Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
end
% >> creates filename and pathname

%% select output folder
% ! have one output folder per experiment
% filenames are not specific to experiment, so important to separate the
% outputs

out_path = uigetdir([], 'Select output folder');
% add last /
out_path = strcat(out_path, '/');

%% import data

% import .mat file
matexp = load(strcat(pathname, filename)); % matexp for .mat file for this experiment

%% write group_tags

csvwrite(strcat(out_path, 'grptags.csv'), matexp.group_tags);

%% write parameter matrix for each day/night
% parameter_matrix is built by Marcus as:
% list of 7 slots, each slot is one day/night
    % 1/ day0
    % 2/ night0
    % 3/ day1
    % 4/ night1
    % 5/ day2
    % 6/ night2
    % 7/ day3
% (day0 and day3 are incomplete so we never use them)
% and each slot is a dataframe rows = fish // columns = parameters
% each datapoint
    % = mean parameter for this day/nigh; this fish; this parameter
% rows: ! only fish included in genotype file are listed
    % match with group_tags for group assignments
% columns: 12 parameters
    % 1/ Active bout length
    % 2/ Active bout mean
    % 3/ Active bout standard deviation
    % 4/ Active bout total
    % 5/ Active bout minimum
    % 6/ Active bout maximum
    % 7/ Number of active bouts
    % 8/ Total time active
    % 9/ Total activity
    % 10/ Inactive bout length
    % 11/ Number of inactive bouts
    % 12/ Total time inactive
% not all are useful/interesting but will export them all here

% to access, command looks like
% parameter_matrix(fish , parameter , day/night)

% write parameter_matrix for each day/night
% starts at 2 and stops at 6 as we skip day0 and day3, see above

% filename are pam_day/night.csv
% pam for parameter_matrix
% day/night will be night0 / day1 / night1 / day2 / night2

csvwrite(strcat(out_path, 'pam_night0.csv'), matexp.parameter_matrix(:, :, 2));
csvwrite(strcat(out_path, 'pam_day1.csv'), matexp.parameter_matrix(:, :, 3));
csvwrite(strcat(out_path, 'pam_night1.csv'), matexp.parameter_matrix(:, :, 4));
csvwrite(strcat(out_path, 'pam_day2.csv'), matexp.parameter_matrix(:, :, 5));
csvwrite(strcat(out_path, 'pam_night2.csv'), matexp.parameter_matrix(:, :, 6));