clear all;
close all;

base_path = 'E:\cmu\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

%% Options Optimized
options_optimized = {};
options_optimized.feature_set = 'e15_flex';
options_optimized.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
options_optimized.use_mirror_motion = false;
options_optimized.use_feature_projection = false;
options_optimized.feature_projection_k = 32;

options_optimized.use_normal_weights = false;
options_optimized.normal_weights_sigma = 4;
options_optimized.k = 800;
options_optimized.generalized_radius = 29.6;
options_optimized.radius = options_optimized.generalized_radius * sqrt( numel(options_optimized.frame_offsets) );
options_optimized.frameRate = 30;

options_optimized.executeChangeFrameRate = true;
options_optimized.executeCreateFeatureSet = true;
options_optimized.executeFindMotionSimilarities = true;
options_optimized.executeFindMotionSegments = true;
options_optimized.executeClusterMotionSegments = false;
options_optimized.generalized_radius_activities = 25;
options_optimized.radius_activities = options_optimized.generalized_radius_activities * sqrt( numel(options_optimized.frame_offsets) );
options_optimized.step_over_frame_penalty = options_optimized.radius;
options_optimized.radius_clustering = options_optimized.radius;


%% Options Optimized Tolerant
options_tolerant = {};
options_tolerant.feature_set = 'e15_flex';
options_tolerant.frame_offsets = [-6 -4 -2 0 2 4 6]; % tolerant best
%options2.frame_offsets = [-15 -10 -5 0 5 10 15]; % strict best
options_tolerant.use_mirror_motion = false;
options_tolerant.use_feature_projection = false;
options_tolerant.feature_projection_k = 32;

options_tolerant.use_normal_weights = false;
options_tolerant.normal_weights_sigma = 4;
options_tolerant.k = 800;
options_tolerant.generalized_radius = 29.6;
options_tolerant.radius = options_tolerant.generalized_radius * sqrt( numel(options_tolerant.frame_offsets) );
options_tolerant.frameRate = 30;

options_tolerant.executeChangeFrameRate = true;
options_tolerant.executeCreateFeatureSet = true;
options_tolerant.executeFindMotionSimilarities = true;
options_tolerant.executeFindMotionSegments = true;
options_tolerant.executeClusterMotionSegments = false;
options_tolerant.generalized_radius_activities = 32;
options_tolerant.radius_activities = options_tolerant.generalized_radius_activities * sqrt( numel(options_tolerant.frame_offsets) );
options_tolerant.step_over_frame_penalty = options_tolerant.radius;
options_tolerant.radius_clustering = options_tolerant.radius;

% % % %% Options Mirrored
% % % options_mirrored = {};
% % % options_mirrored.feature_set = 'e15_flex';
% % % options_mirrored.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
% % % options_mirrored.use_mirror_motion = true;
% % % options_mirrored.use_feature_projection = false;
% % % options_mirrored.feature_projection_k = 32;
% % % 
% % % options_mirrored.use_normal_weights = false;
% % % options_mirrored.normal_weights_sigma = 4;
% % % options_mirrored.k = 800;
% % % options_mirrored.generalized_radius = 29.6;
% % % options_mirrored.radius = options_mirrored.generalized_radius * sqrt( numel(options_mirrored.frame_offsets) );
% % % options_mirrored.frameRate = 30;
% % % 
% % % options_mirrored.executeChangeFrameRate = true;
% % % options_mirrored.executeCreateFeatureSet = true;
% % % options_mirrored.executeFindMotionSimilarities = true;
% % % options_mirrored.executeFindMotionSegments = true;
% % % options_mirrored.executeClusterMotionSegments = false;
% % % options_mirrored.generalized_radius_activities = 25;
% % % options_mirrored.radius_activities = options_mirrored.generalized_radius_activities * sqrt( numel(options_mirrored.frame_offsets) );
% % % options_mirrored.step_over_frame_penalty = options_mirrored.radius;
% % % options_mirrored.radius_clustering = options_mirrored.radius;


%% Options Mirrored
options_mirrored = {};
options_mirrored.feature_set = 'e15_flex';
options_mirrored.frame_offsets = [-10 -5 0 5 10];
options_mirrored.use_mirror_motion = true;
options_mirrored.use_feature_projection = false;
options_mirrored.feature_projection_k = 32;

options_mirrored.use_normal_weights = false;
options_mirrored.normal_weights_sigma = 4;
options_mirrored.k = 800;
options_mirrored.generalized_radius = 30;
options_mirrored.radius = options_mirrored.generalized_radius * sqrt( numel(options_mirrored.frame_offsets) );
options_mirrored.frameRate = 30;

options_mirrored.executeChangeFrameRate = true;
options_mirrored.executeCreateFeatureSet = true;
options_mirrored.executeFindMotionSimilarities = true;
options_mirrored.executeFindMotionSegments = true;
options_mirrored.executeClusterMotionSegments = false;
options_mirrored.generalized_radius_activities = 25;
options_mirrored.radius_activities = options_mirrored.generalized_radius_activities * sqrt( numel(options_mirrored.frame_offsets) );
options_mirrored.step_over_frame_penalty = options_mirrored.radius;
options_mirrored.radius_clustering = options_mirrored.radius;

%% Options Bundled Features
options_bundled = {};
options_bundled.feature_set = 'e15_flex';
options_bundled.frame_offsets = [-10 -5 0 5 10];
options_bundled.use_mirror_motion = false;
options_bundled.use_feature_projection = true;
options_bundled.feature_projection_k = 32;

options_bundled.use_normal_weights = false;
options_bundled.normal_weights_sigma = 4;
options_bundled.k = 800;
options_bundled.generalized_radius = 21;
options_bundled.radius = options_bundled.generalized_radius * sqrt( numel(options_bundled.frame_offsets) );
options_bundled.frameRate = 30;

options_bundled.executeChangeFrameRate = true;
options_bundled.executeCreateFeatureSet = true;
options_bundled.executeFindMotionSimilarities = true;
options_bundled.executeFindMotionSegments = true;
options_bundled.executeClusterMotionSegments = false;
options_bundled.generalized_radius_activities = 11;
options_bundled.radius_activities = options_bundled.generalized_radius_activities * sqrt( numel(options_bundled.frame_offsets) );
options_bundled.step_over_frame_penalty = options_bundled.radius;
options_bundled.radius_clustering = options_bundled.radius;

%% Options Mirrored Bundled Features
options_bundled_mirrored = {};
options_bundled_mirrored.feature_set = 'e15_flex';
options_bundled_mirrored.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
options_bundled_mirrored.use_mirror_motion = true;
options_bundled_mirrored.use_feature_projection = true;
options_bundled_mirrored.feature_projection_k = 32;

options_bundled_mirrored.use_normal_weights = false;
options_bundled_mirrored.normal_weights_sigma = 4;
options_bundled_mirrored.k = 800;
options_bundled_mirrored.generalized_radius = 25;
options_bundled_mirrored.radius = options_bundled_mirrored.generalized_radius * sqrt( numel(options_bundled_mirrored.frame_offsets) );
options_bundled_mirrored.frameRate = 30;

options_bundled_mirrored.executeChangeFrameRate = true;
options_bundled_mirrored.executeCreateFeatureSet = true;
options_bundled_mirrored.executeFindMotionSimilarities = true;
options_bundled_mirrored.executeFindMotionSegments = true;
options_bundled_mirrored.executeClusterMotionSegments = false;
options_bundled_mirrored.generalized_radius_activities = 11;
options_bundled_mirrored.radius_activities = options_bundled_mirrored.generalized_radius_activities * sqrt( numel(options_bundled_mirrored.frame_offsets) );
options_bundled_mirrored.step_over_frame_penalty = options_bundled_mirrored.radius;
options_bundled_mirrored.radius_clustering = options_bundled_mirrored.radius;


%options_list = {options_optimized, options_tolerant, options_mirrored};
%record_names = {'Optimized', 'Optimized Tolerant', 'Optimized Mirrored'};

options_list = {options_optimized, options_tolerant, options_bundled, options_mirrored, options_bundled_mirrored};
record_names = {'Optimized', 'Optimized Tolerant', 'Bundled Features', 'Optimized Mirrored', 'Bundled + Mirrored'};

range = 1:14;
%range = 1:2;

range_count = numel(range);
options_list_count = numel(options_list);
names = cell(1, range_count);
frame_counts = zeros(range_count, options_list_count);
strict_results = zeros(range_count, options_list_count);
tolerant_results = zeros(range_count, options_list_count);

[ gt, annotation ] = tw_gt_86();
subcutsHACA = tw_get_subcuts_haca();
subcutsACA = tw_get_subcuts_aca();
clusteringZhouHACA = tw_get_clustering_haca();
clusteringZhouACA = tw_get_clustering_aca();

for trial=range
    gt_data = gt{trial};
    annotation_data = annotation{trial};
    [ gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color ] = tw_gt_to_segments(gt_data, annotation_data);
    [ skel, mot, name ] = tw_get_mot_by_index(trial, base_path, skel_file, file_prefix, file_suffix);
    
    % + gt + old + HACA + ACA
    segment_list_count = options_list_count + 4;
    
    segments_start_list = cell(1, segment_list_count);
    segments_end_list = cell(1, segment_list_count);
    segments_color_list = cell(1, segment_list_count);
    segments_cuts_list = cell(1, segment_list_count);
    segments_subcuts_list = cell(1, segment_list_count);
    
    segments_start_list{1} = gt_segment_start;
    segments_end_list{1} = gt_segment_end;
    segments_color_list{1} = gt_segment_color;
    
%     frame_count = gt_segment_end(end);
    
    % Old Implementation
    mot_old = changeFrameRate(skel, mot, 30);
    [mots, submots, comps, cuts, allsubcuts] = tw_LNGConnectedComponents(skel, mot_old);
    frame_count = mot_old.nframes;
    [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
    [ strict, tolerant, cut_segment_to_gt_segment, cut_segment_color ] = tw_segment_match( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, cut_segment_start, cut_segment_end);

    segments_start_list{end - 2} = cut_segment_start;
    segments_end_list{end - 2} = cut_segment_end;
    segments_color_list{end - 2} = cut_segment_color;
    segments_cuts_list{end - 2} = cuts;
    segments_subcuts_list{end - 2} = allsubcuts;
    
    % HACA
    [ haca_segment_start, haca_segment_end ] = tw_cuts_to_segments( subcutsHACA{trial}, frame_count );
    [ haca_segment_color ] = tw_getZhouSegmentColors( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, haca_segment_start, haca_segment_end, clusteringZhouHACA{trial} );
    segments_start_list{end - 1} = haca_segment_start;
    segments_end_list{end - 1} = haca_segment_end;
    segments_color_list{end - 1} = haca_segment_color;
    
    % ACA
    [ aca_segment_start, aca_segment_end ] = tw_cuts_to_segments( subcutsACA{trial}, frame_count );
    [ aca_segment_color ] = tw_getZhouSegmentColors( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, aca_segment_start, aca_segment_end, clusteringZhouACA{trial} );
    segments_start_list{end} = aca_segment_start;
    segments_end_list{end} = aca_segment_end;
    segments_color_list{end} = aca_segment_color;
    
    for options_list_idx=1:options_list_count
        options_optimized = options_list{options_list_idx};

        [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options_optimized);
        frame_count = mot.nframes;
        [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
        [ strict, tolerant, cut_segment_to_gt_segment, cut_segment_color ] = tw_segment_match( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, cut_segment_start, cut_segment_end);

        names{trial} = name;
        frame_counts(trial, options_list_idx) = frame_count;
        strict_results(trial, options_list_idx) = strict;
        tolerant_results(trial, options_list_idx) = tolerant;
        
        segments_start_list{options_list_idx + 1} = cut_segment_start;
        segments_end_list{options_list_idx + 1} = cut_segment_end;
        segments_color_list{options_list_idx + 1} = cut_segment_color;
        segments_cuts_list{options_list_idx + 1} = cuts;
        segments_subcuts_list{options_list_idx + 1} = subcuts;
    end
    
    figure(trial);
    compare_names = cell(1, numel(record_names) + 4);
    compare_names{1} = 'Ground Truth';    
    compare_names(2:(end - 3)) = record_names;
    compare_names{end - 2} = 'Vögele et al. 2014';
    compare_names{end - 1} = 'HACA';
    compare_names{end} = 'ACA';
    tw_plotCompareSegmentation(segments_cuts_list, segments_subcuts_list, segments_start_list, segments_end_list, segments_color_list, compare_names, names{trial});   
end

[ strict_old, tolerant_old, old_record_names ] = tw_old_evaluation();

record_names_ext = cell(1, numel(record_names) + numel(old_record_names));
for j = 1:numel(record_names)
    record_names_ext{j} = record_names{j};
end
for j = 1:numel(old_record_names)
    record_names_ext{numel(record_names) + j} = old_record_names{j};
end
figure(101);
strict_old_range = strict_old(range, :);
strict_records_count = size(strict_results, 2);
strict_old_records_count = size(strict_old_range, 2);
strict_ext = zeros(range_count, strict_records_count + strict_old_records_count);
strict_ext(:, 1:strict_records_count) = strict_results;  
strict_ext(:, (strict_records_count + 1):end) = strict_old_range;
strict_all = tw_accuracyPlot( strict_ext, 'Strict Evaluation', names, record_names_ext, tw_segment_colors(), frame_counts );
figure(102);
tolerant_old_range = tolerant_old(range, :);
tolerant_records_count = size(tolerant_results, 2);
tolerant_old_records_count = size(tolerant_old_range, 2);
tolerant_ext = zeros(range_count, tolerant_records_count + tolerant_old_records_count);
tolerant_ext(:, 1:tolerant_records_count) = tolerant_results;  
tolerant_ext(:, (tolerant_records_count + 1):end) = tolerant_old_range;
tolerant_all = tw_accuracyPlot( tolerant_ext, 'Tolerant Evaluation', names, record_names_ext, tw_segment_colors(), frame_counts );
