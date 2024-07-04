function [comps, sframes, eframes] = call_segmentation(tag)
%CALL_SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

str = num2str(tag);
    
% Determine the filename based on the value of x
if tag < 10
    filename = strcat('86_0', str, '.amc');
else
    filename = strcat('86_', str, '.amc');
end

[skel, mot] = readMocap(['../CMU_86/','86.asf'], ['../CMU_86/',filename], [], true, true, true, 0);
% [ skel, mot ] = readMocapGUI();

% go back to the source code.

% set up a simple options struct.
options.DEBUG = false;

% some more debugging levels:
options.DEBUG_FIND_MOTION_ACTIVITIES = false;
options.DEBUG_FIND_MOTION_SEGMENTS = false;
options.DEBUG_CLUSTER_MOTION_SEGMENTS = false;

% all other options set in tw_segmentation and the segmentation methods can
% be set here, and will override the default options given in the code.

%option.MYOPTION = XXX;

% start segmentation
[mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);

n = numel(submots);
sframes = zeros(n, 1);
eframes = zeros(n, 1);

% Loop through each cell to extract the nframes property
for i = 1:n
    % Check if the current cell contains a struct
    if ~isstruct(submots{i})
        error('submots{%d} is not a struct', i);
    end

    % Check if the struct contains the 'sf' property
    if ~isfield(submots{i}, 'sf')
        error('The struct in submots{%d} does not contain the sf property', i);
    end

    % Extract the ef property
    sframes(i) = submots{i}.sf;

    if ~isfield(submots{i}, 'ef')
        error('The struct in submots{%d} does not contain the ef property', i);
    end

    % Extract the ef property
    eframes(i) = submots{i}.ef;
end

end

