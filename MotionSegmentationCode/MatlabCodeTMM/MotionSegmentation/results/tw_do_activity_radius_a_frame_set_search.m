clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

frame_offset_list_s5 = {
    'e15_flex', 0;
    'e15_flex', [-5 0 5];
    'e15_flex', [-10 -5 0 5 10];
    'e15_flex', [-15 -10 -5 0 5 10 15];
    'e15_flex', [-20 -15 -10 -5 0 5 10 15 20];
    'e15_flex', [-25 -20 -15 -10 -5 0 5 10 15 20 25]
};

frame_offset_list_e30_s5 = {
    'e30_flex', 0;
    'e30_flex', [-5 0 5];
    'e30_flex', [-10 -5 0 5 10];
    'e30_flex', [-15 -10 -5 0 5 10 15]
};

frame_offset_list_s5u = {
    'e15_flex', [-5 0];
    'e15_flex', [-10 -5 0];
    'e15_flex', [-15 -10 -5 0];
    'e15_flex', [0 5];
    'e15_flex', [0 5 10];
    'e15_flex', [0 5 10 15]
};

frame_offset_list_s3 = {
    'e15_flex', [-3 0 3];
    'e15_flex', [-6 -3 0 3 6];
    'e15_flex', [-9 -6 -3 0 3 6 9];
    'e15_flex', [-12 -9 -6 -3 0 3 6 9 12];
    'e15_flex', [-15 -12 -9 -6 -3 0 3 6 9 12 15]
};

frame_offset_list_s2 = {
    'e15_flex', [-2 0 2];
    'e15_flex', [-4 -2 0 2 4];
    'e15_flex', [-6 -4 -2 0 2 4 6];
    'e15_flex', [-8 -6 -4 -2 0 2 4 6 8];
    'e15_flex', [-10 -8 -6 -4 -2 0 2 4 6 8 10]
};

frame_offset_list_e30_s2 = {
    'e30_flex', [-2 0 2];
    'e30_flex', [-4 -2 0 2 4];
    'e30_flex', [-6 -4 -2 0 2 4 6];
    'e30_flex', [-8 -6 -4 -2 0 2 4 6 8];
    'e30_flex', [-10 -8 -6 -4 -2 0 2 4 6 8 10]
};

frame_offset_list_s2_8 = {
    'e15_flex', [-8 0 8];
    'e15_flex', [-8 -4 0 4 8];
    'e15_flex', [-8 -6 -4 -2 0 2 4 6 8];
    'e15_flex', [-8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
};

frame_offset_list_best_a = {
    'e15_flex', [-5 0 5];
    'e15_flex', [-10 -5 0 5 10];
    'e15_flex', [-9 -6 -3 0 3 6 9];
    'e15_flex', [-8 -6 -4 -2 0 2 4 6 8];
    'e15_flex', [-15 -10 -5 0];
    'e15_flex', [0 5 10 15]
};

frame_offset_list_e30_best_a = {
    'e30_flex', [-5 0 5];
    'e30_flex', [-10 -5 0 5 10];
    'e30_flex', [-9 -6 -3 0 3 6 9];
    'e30_flex', [-8 -6 -4 -2 0 2 4 6 8];
    'e30_flex', [-15 -10 -5 0];
    'e30_flex', [0 5 10 15]
};

frame_offset_list = frame_offset_list_best_a;

is_e30 = false;
use_feature_projection = true;

frame_offset_count = size(frame_offset_list, 1);

trial_range = 1:14;
%trial_range = 1:1;

% never filter through k
%k = 1184; % radius = 42
%k = 1039; % radius = 34 bundled features
k = 1233; % radius = 40 bundled features
if (use_feature_projection)
    ref_radius_a = 22;
    offset_radius_a = 18;
    step_radius_a = 1;
else
    ref_radius_a = 30;
    offset_radius_a = 12;
    step_radius_a = 1;
end
radius_a_range = (ref_radius_a - offset_radius_a):step_radius_a:(ref_radius_a + offset_radius_a);
radius = max(radius_a_range);

[ gt, annotation ] = tw_gt_86();

strict_results = zeros(frame_offset_count, numel(trial_range), numel(radius_a_range));
tolerant_results = zeros(frame_offset_count, numel(trial_range), numel(radius_a_range));
activity_count_results = zeros(frame_offset_count, numel(trial_range), numel(radius_a_range));

for feature_idx=1:frame_offset_count
    feature_set = frame_offset_list{feature_idx, 1};
    frame_offsets = frame_offset_list{feature_idx, 2};
    
    for trial_idx=1:numel(trial_range)
        trial = trial_range(trial_idx);
        [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );

        gt_trial = gt{trial};
        annotation_trial = annotation{trial};
        [ gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color ] = tw_gt_to_segments(gt_trial, annotation_trial);

        pre_options.feature_set = feature_set;
        pre_options.frame_offsets = frame_offsets;
        pre_options.use_mirror_motion = false;
        pre_options.use_feature_projection = use_feature_projection;
        pre_options.feature_projection_k = 32;
        pre_options.use_normal_weights = false;
        pre_options.normal_weights_sigma = 4;
        pre_options.k = k;
        pre_options.generalized_radius = radius;
        pre_options.radius = pre_options.generalized_radius * sqrt( numel(pre_options.frame_offsets) );
        pre_options.frameRate = 30;

        mot = changeFrameRate(skel, mot, pre_options.frameRate);
        pre_options.k = min(mot.nframes, pre_options.k);

        [ fmat, fmat_mirror ] = tw_createFeatureSet( skel, mot, pre_options );
        [ nnidx, nndists, nnidx_mirror, nndists_mirror ] = tw_findMotionSimilarities( fmat, fmat_mirror, pre_options );

        for radius_idx=1:numel(radius_a_range)
            radius_a = radius_a_range(radius_idx);

            options = {};
            options.frameRate = pre_options.frameRate;
            options.feature_set = pre_options.feature_set;
            options.frame_offsets = pre_options.frame_offsets;
            options.use_mirror_motion = pre_options.use_mirror_motion;
            options.use_feature_projection = pre_options.use_feature_projection;
            options.feature_projection_k = pre_options.feature_projection_k;
            options.use_normal_weights = pre_options.use_normal_weights;
            options.normal_weights_sigma = pre_options.normal_weights_sigma;
            options.k = pre_options.k;
            options.generalized_radius = pre_options.generalized_radius;
            options.radius = pre_options.radius;

            options.executeChangeFrameRate = false;
            options.executeCreateFeatureSet = false;
            options.executeFindMotionSimilarities = false;
            options.nnidx = nnidx;
            options.nndists = nndists;
            options.nnidx_mirror = nnidx_mirror;
            options.nndists_mirror = nndists_mirror;
            options.executeFindMotionSegments = false;
            options.executeClusterMotionSegments = false;
            options.generalized_radius_activities = radius_a;
            options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
            options.step_over_frame_penalty = options.radius;
            options.radius_clustering = options.radius;

            [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);

            frame_count = mot.nframes;
            [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
            [ strict, tolerant, cut_segment_to_gt_segment, cut_segment_color ] = tw_segment_match( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, cut_segment_start, cut_segment_end);

            strict_results(feature_idx, trial_idx, radius_idx) = strict;
            tolerant_results(feature_idx, trial_idx, radius_idx) = tolerant;
            activity_count_results(feature_idx, trial_idx, radius_idx) = numel(cut_segment_start) / 20;
        end
    end
end

strict_eval = squeeze(sum(strict_results, 2)) / numel(trial_range);
tolerant_eval = squeeze(sum(tolerant_results, 2)) / numel(trial_range);
activity_eval = squeeze(sum(activity_count_results, 2)) / numel(trial_range);


%% visualization
feature_set_names = cell(1, frame_offset_count);
for i=1:frame_offset_count
    feature_set_names{i} = strcat('[',  strrep(strtrim(sprintf('%d ', frame_offset_list{i, 2})), ' ', ','), ']');
end

contour_range = 0:0.01:0.99;
caxis_value = [0.74 1];
colorbar_ticks = caxis_value(1):0.02:1;
feature_range = 1:frame_offset_count;

figure(1);
ax = gca;
[C,h] = contourf(radius_a_range, feature_range, strict_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Feature Set Offsets');
[void, max_s_index] = max(strict_eval(:));
[max_s_row, max_s_column] = ind2sub(size(strict_eval), max_s_index);
max_s_radius = radius_a_range(max_s_column);
if (use_feature_projection)
    title('Strict Accuracy Bundled F^{15}_E Feature Sets');
else
    if (is_e30)
        title('Strict Accuracy F^{30}_E Feature Sets');
    else
        title('Strict Accuracy F^{15}_E Feature Sets');
    end
end

figure(2);
ax = gca;
[C,h] = contourf(radius_a_range, feature_range, tolerant_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Feature Set Offsets');
[void, max_t_index] = max(tolerant_eval(:));
[max_t_row, max_t_column] = ind2sub(size(tolerant_eval), max_t_index);
max_t_radius = radius_a_range(max_t_column);
if (use_feature_projection)
    title('Tolerant Accuracy Bundled F^{15}_E Feature Sets');
else
    if (is_e30)
        title('Tolerant Accuracy F^{30}_E Feature Sets');
    else
        title('Tolerant Accuracy F^{15}_E Feature Sets');
    end
end

figure(3);
ax = gca;
[C,h] = contourf(radius_a_range, feature_range, activity_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Feature Set Offsets');
[void, max_a_index] = max(activity_eval(:));
[max_a_row, max_a_column] = ind2sub(size(activity_eval), max_a_index);
max_a_radius = radius_a_range(max_a_column);
if (use_feature_projection)
    title('Relative Activity Count Bundled F^{15}_E Feature Sets');  
else
    if (is_e30)
        title('Relative Activity Count F^{30}_E Feature Sets');
    else
        title('Relative Activity Count F^{15}_E Feature Sets');
    end
end

fprintf('Strict Maximum := (Radius: %f, Feature Set: %s, Value: %f)\n', max_s_radius, feature_set_names{max_s_row}, strict_eval(max_s_row, max_s_column));
fprintf('Tolerant Maximum := (Radius: %f, Feature Set: %s, Value: %f)\n', max_t_radius, feature_set_names{max_t_row}, tolerant_eval(max_t_row, max_t_column));
fprintf('Relative Activity Count Maximum := (Radius: %f, Feature Set: %s, Value: %f)\n', max_a_radius, feature_set_names{max_a_row}, activity_eval(max_a_row, max_a_column));