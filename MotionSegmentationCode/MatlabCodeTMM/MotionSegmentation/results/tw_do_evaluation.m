clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

%% Options 1
options = {};
options.feature_set = 'e15_flex';
options.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
options.use_mirror_motion = false;
options.use_feature_projection = false;
options.feature_projection_k = 32;

options.use_normal_weights = false;
options.normal_weights_sigma = 4;
options.k = 800;
options.generalized_radius = 29.6;
options.radius = options.generalized_radius * sqrt( numel(options.frame_offsets) );
options.frameRate = 30;

options.executeChangeFrameRate = true;
options.executeCreateFeatureSet = true;
options.executeFindMotionSimilarities = true;
options.executeFindMotionSegments = false;
options.executeClusterMotionSegments = false;
options.generalized_radius_activities = 25;
options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
options.step_over_frame_penalty = options.radius;
options.radius_clustering = options.radius;


%% Options 2
options2 = {};
options2.feature_set = 'e15_flex';
options2.frame_offsets = [-6 -4 -2 0 2 4 6]; % tolerant best
%options2.frame_offsets = [-15 -10 -5 0 5 10 15]; % strict best
options2.use_mirror_motion = false;
options2.use_feature_projection = false;
options2.feature_projection_k = 32;

options2.use_normal_weights = false;
options2.normal_weights_sigma = 4;
options2.k = 800;
options2.generalized_radius = 29.6;
options2.radius = options2.generalized_radius * sqrt( numel(options2.frame_offsets) );
options2.frameRate = 30;

options2.executeChangeFrameRate = true;
options2.executeCreateFeatureSet = true;
options2.executeFindMotionSimilarities = true;
options2.executeFindMotionSegments = false;
options2.executeClusterMotionSegments = false;
options2.generalized_radius_activities = 32;
options2.radius_activities = options2.generalized_radius_activities * sqrt( numel(options2.frame_offsets) );
options2.step_over_frame_penalty = options2.radius;
options2.radius_clustering = options2.radius;


%options_list = {options1};
options_list = {options, options2};
%record_names = {'Our Algorithm'};
record_names = {'Optimized', 'Optimized Tolerant'};


range = 1:14;
%range = 1:2;
is_parameter_search = false;
tw_evaluation(range, options_list, record_names, base_path, skel_file, file_prefix, file_suffix);