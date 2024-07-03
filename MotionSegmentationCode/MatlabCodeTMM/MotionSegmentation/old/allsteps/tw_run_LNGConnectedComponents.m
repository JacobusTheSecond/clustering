close all;
clear all;

% standard
cd 'D:\mocap\cmu\all_asfamc\subjects\86';
% dance
%cd 'D:\mocap\cmu\all_asfamc\subjects\05';
% vignettes - locomotion, upper-body motions (focus: motion transitions)
%cd 'D:\mocap\cmu\all_asfamc\subjects\56';
% modern dance, gymnastics
%cd 'D:\mocap\cmu\all_asfamc\subjects\49';

[ skel, mot ] = readMocapGUI();
cd( fullfile('D:\SVN\MotionToolsMatlab\', 'theses', 'MotionSegmentation') )
%[skel,mot] = readMocap('D:\mocap\hdm\HDM_bd.asf','D:\mocap\hdm\HDM_bd_01-01_01_120.amc');

% Bewegung resamplen auf 30fps
mot = changeFrameRate(skel,mot,30);

options.DEBUG = true;

%Segmentierung starten
[mots, submots, comps, cuts, allsubcuts] = tw_LNGConnectedComponents(skel, mot, options);