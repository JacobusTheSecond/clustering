clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

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
use_mirror_motion = true;
use_feature_projection = true;

frame_offset_count = size(frame_offset_list, 1);

trial_range = 1:14;
%trial_range = 1:1;

frame_rate = 30;
% never filter through k
%k = 1184; % radius = 42
%k = 949; % radius = 36
%k = 992; % radius = 37
%k = 687; % radius = 25 bundled
k = 1094; % radius = 36 bundled

%k = 965; % radius = 41; fe_30

%radius_a = 24;
%radius_a = 25;
%radius_a = 29;
radius_a = 11;

ref_radius = 24;%radius_a;
offset_radius = 12;%14;%12;
step_radius = 1;
radius_range = ref_radius:step_radius:(ref_radius + offset_radius);


deviation_results = zeros(frame_offset_count, numel(trial_range), numel(radius_range));
segmentation_ratio_results = zeros(frame_offset_count, numel(trial_range), numel(radius_range));
segment_count_results = zeros(frame_offset_count, numel(trial_range), numel(radius_range));

for feature_idx=1:frame_offset_count
    feature_set = frame_offset_list{feature_idx, 1};
    frame_offsets = frame_offset_list{feature_idx, 2};
    
    for trial_idx=1:numel(trial_range)
        trial = trial_range(trial_idx);
        [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );

        pre_options.feature_set = feature_set;
        pre_options.frame_offsets = frame_offsets;
        pre_options.use_mirror_motion = use_mirror_motion;
        pre_options.use_feature_projection = use_feature_projection;
        pre_options.feature_projection_k = 32;

        mot = changeFrameRate(skel, mot, frame_rate);
        k = min(mot.nframes, k);

        [ fmat, fmat_mirror ] = tw_createFeatureSet( skel, mot, pre_options );

        for radius_idx=1:numel(radius_range)
            radius = radius_range(radius_idx);

            options = {};
            options.feature_set = pre_options.feature_set;
            options.frame_offsets = pre_options.frame_offsets;
            options.use_mirror_motion = pre_options.use_mirror_motion;
            options.use_feature_projection = pre_options.use_feature_projection;
            options.feature_projection_k = pre_options.feature_projection_k;
            
            options.use_normal_weights = false;
            options.normal_weights_sigma = 4;
            options.k = k;
            options.generalized_radius = radius;
            options.radius = options.generalized_radius * sqrt( numel(options.frame_offsets) );
            options.frameRate = frame_rate;
            
            options.executeChangeFrameRate = false;
            options.executeCreateFeatureSet = false;
            options.fmat = fmat;
            options.fmat_mirror = fmat_mirror;
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

            deviation_results(feature_idx, trial_idx, radius_idx) = deviation;
            segmentation_ratio_results(feature_idx, trial_idx, radius_idx) = segmentation_ratio;
            segment_count_results(feature_idx, trial_idx, radius_idx) = segment_count / 75;
        end
    end
end

deviation_eval = squeeze(sum(deviation_results, 2)) / numel(trial_range);
segmentation_ratio_eval = squeeze(sum(segmentation_ratio_results, 2)) / numel(trial_range);
segment_count_eval = squeeze(sum(segment_count_results, 2)) / numel(trial_range);

save('backup');

%% visualization
feature_set_names = cell(1, frame_offset_count);
for i=1:frame_offset_count
    feature_set_names{i} = strcat('[',  strrep(strtrim(sprintf('%d ', frame_offset_list{i, 2})), ' ', ','), ']');
end

contour_range = 0:0.01:0.99;
feature_range = 1:frame_offset_count;

figure(1);
ax = gca;
if (use_mirror_motion)
    caxis_value = [0.5 1];
else
    caxis_value = [0 0.5];
end
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
[C,h] = contourf(radius_range, feature_range, deviation_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius');
ylabel('Feature Set Offsets');
[void, min_d_index] = min(deviation_eval(:));
[min_d_row, min_d_column] = ind2sub(size(deviation_eval), min_d_index);
min_d_radius = radius_range(min_d_column);
if (use_feature_projection)
    title('Relative Deviation Bundled F^{15}_E Feature Sets');
else
    if (is_e30)
        title('Relative Deviation F^{30}_E Feature Sets');
    else
        title('Relative Deviation F^{15}_E Feature Sets');
    end
end

figure(2);
ax = gca;
caxis_value = [0.5 1];
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
[C,h] = contourf(radius_range, feature_range, segmentation_ratio_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius');
ylabel('Feature Set Offsets');
[void, max_s_index] = max(segmentation_ratio_eval(:));
[max_s_row, max_s_column] = ind2sub(size(segmentation_ratio_eval), max_s_index);
max_s_radius = radius_range(max_s_column);
if (use_feature_projection)
    title('Segmentation Ratio Bundled F^{15}_E Feature Sets');
else
    if (is_e30)
        title('Segmentation Ratio F^{30}_E Feature Sets');
    else
        title('Segmentation Ratio F^{15}_E Feature Sets');
    end
end

figure(3);
ax = gca;
if (use_mirror_motion)
    caxis_value = [0.5 1];
else
    caxis_value = [0.5 0.75];
end
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
[C,h] = contourf(radius_range, feature_range, segment_count_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_range);
set(ax, 'YTick', feature_range);
set(ax, 'YTickLabel', feature_set_names);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius');
ylabel('Feature Set Offsets');
[void, max_r_index] = max(segment_count_eval(:));
[max_r_row, max_r_column] = ind2sub(size(segment_count_eval), max_r_index);
max_r_radius = radius_range(max_r_column);
if (use_feature_projection)
    title('Relative Segment Count Bundled F^{15}_E Feature Sets');
else
    if (is_e30)
        title('Relative Segment Count F^{30}_E Feature Sets');
    else
        title('Relative Segment Count F^{15}_E Feature Sets');
    end
end

fprintf('Deviation Minimum := (Radius: %f, Feature Set: %s, Value: %f)\n', min_d_radius, feature_set_names{min_d_row}, deviation_eval(min_d_row, min_d_column));
fprintf('Segmentation Ratio Maximum := (Radius: %f, Feature Set: %s, Value: %f)\n', max_s_radius, feature_set_names{max_s_row}, segmentation_ratio_eval(max_s_row, max_s_column));
fprintf('Relative Segment Count := (Radius: %f, Feature Set: %s, Value: %f)\n', max_r_radius, feature_set_names{max_r_row}, segment_count_eval(max_r_row, max_r_column));