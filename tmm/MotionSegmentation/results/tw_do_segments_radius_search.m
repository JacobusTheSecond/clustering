clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

feature_set = 'e15_flex';
frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
use_mirror_motion = true;

trial_range = 1:14;
%trial_range = 12:12;
%trial_range = 4:4;

k = 755; % radius = 31
radius_a = 25;
ref_radius = 29;
offset_radius = 2;
step_radius = 0.1;
radius_range = ref_radius:step_radius:(ref_radius + offset_radius);

title_name = strcat('Radius Activity:{ }', num2str(radius_a), '{ | }Feature Set: F^{15}_E [',  strrep(strtrim(sprintf('%d ', frame_offsets)), ' ', ','), ']');

deviation_results = zeros(numel(trial_range), numel(radius_range));
ratio_results = zeros(numel(trial_range), numel(radius_range));
segment_count_results = zeros(numel(trial_range), numel(radius_range));


for radius_idx=1:numel(radius_range)
    radius = radius_range(radius_idx);
    
    for trial_idx=1:numel(trial_range)
        trial = trial_range(trial_idx);
        [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );
        
        options = {};
        options.feature_set = feature_set;
        options.frame_offsets = frame_offsets;
        options.use_mirror_motion = use_mirror_motion;
        options.use_feature_projection = false;
        options.feature_projection_k = 32;
        options.use_normal_weights = false;
        options.normal_weights_sigma = 4;
        options.k = k;
        options.generalized_radius = radius;
        options.radius = options.generalized_radius * sqrt( numel(options.frame_offsets) );
        options.frameRate = 30;

        mot = changeFrameRate(skel, mot, options.frameRate);
        options.k = min(mot.nframes, options.k);
 
        options.executeChangeFrameRate = false;
        options.executeCreateFeatureSet = true;
        options.executeFindMotionSimilarities = true;
        options.executeFindMotionSegments = true;
        options.executeClusterMotionSegments = false;
        options.generalized_radius_activities = radius_a;
        options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
        options.step_over_frame_penalty = options.radius;
        options.radius_clustering = options.radius;

        [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);

        frame_count = mot.nframes;
        [ deviation, segmentation_ratio, activity_count, segment_count ] = tw_relative_segment_deviation(cuts, subcuts, frame_count);
        
        deviation_results(trial_idx, radius_idx) = deviation;
        ratio_results(trial_idx, radius_idx) = segmentation_ratio;
        segment_count_results(trial_idx, radius_idx) = segment_count / 75;
    end
end

deviation_eval = sum(deviation_results, 1) / numel(trial_range);
ratio_eval = sum(ratio_results, 1) / numel(trial_range);
segment_count_eval = sum(segment_count_results, 1) / numel(trial_range);

%% visualization

figure(1);
plot(radius_range, deviation_eval, 'r', radius_range, segment_count_eval, 'b', 'LineWidth', 2);
xlabel('Radius');
legend('Relative Deviation', 'Relative Segment Count');
ax = gca;
set(ax, 'XTick', radius_range);
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'on');
title(title_name);

figure(2);
plot(radius_range, deviation_eval, 'r', 'LineWidth', 2);
xlabel('Radius');
ylabel('Relative Deviation')
ax = gca;
set(ax, 'XTick', radius_range);
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'on');
title(title_name);