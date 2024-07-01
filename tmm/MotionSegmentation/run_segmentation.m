close all;
clear all;

% move to the directory with the example mocap data.
cd '..\CMU_86';


% use file reader to open asf/amc file
[ skel, mot ] = readMocapGUI();

% go back to the source code.
cd( '..\MotionSegmentation' )

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