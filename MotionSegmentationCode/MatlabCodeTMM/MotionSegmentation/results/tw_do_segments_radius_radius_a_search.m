clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

feature_set = 'e15_flex';
frame_offsets = [-5 0 5];
% 0;
% [-5 0 5];
% [-10 -5 0 5 10];
feature_set_name = strcat('[',  strrep(strtrim(sprintf('%d ', frame_offsets)), ' ', ','), ']');
feature_set_name = strcat('(Feature Set F^{15}_E{ }', feature_set_name, ')');
frameRate = 30;

%trial_range = 1:14;
trial_range = 1:1;

% never filter through k
%k = 931; % for radius = 38
%k = 1101; % for radius = 41
k = 695; % for radius = 34

ref_radius_a = 20;
offset_radius_a = 10;
step_radius_a = 1;
radius_a_range = ref_radius_a:step_radius_a:(ref_radius_a + offset_radius_a);

ref_radius = ref_radius_a;
offset_radius = 14;
step_radius = 1;
radius_range = ref_radius:step_radius:(ref_radius + offset_radius);

deviation_results = zeros(numel(trial_range), numel(radius_range), numel(radius_a_range));
segmentation_ratio_results = zeros(numel(trial_range), numel(radius_range), numel(radius_a_range));
activity_count_results = zeros(numel(trial_range), numel(radius_range), numel(radius_a_range));

for trial_idx=1:numel(trial_range)
    trial = trial_range(trial_idx);
    [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );

    mot = changeFrameRate(skel, mot, frameRate);
    k = min(mot.nframes, k);
    
    pre_options.feature_set = feature_set;
    pre_options.frame_offsets = frame_offsets;
    pre_options.use_mirror_motion = false;
    pre_options.use_feature_projection = false;
    pre_options.feature_projection_k = 32;
    pre_options.frameRate = frameRate;
    
    [ fmat, fmat_mirror ] = tw_createFeatureSet( skel, mot, pre_options );
    
    for radius_idx = 1:numel(radius_range)
        radius = radius_range(radius_idx);
        
        pre_options.frame_offsets = frame_offsets;
        pre_options.use_normal_weights = false;
        pre_options.normal_weights_sigma = 4;
        pre_options.k = k;
        pre_options.generalized_radius = radius;
        pre_options.radius = pre_options.generalized_radius * sqrt( numel(pre_options.frame_offsets) );

        [ nnidx, nndists, nnidx_mirror, nndists_mirror ] = tw_findMotionSimilarities( fmat, fmat_mirror, pre_options );

        for radius_a_idx=1:numel(radius_a_range)
            radius_a = radius_a_range(radius_a_idx);
            
            if (radius_a > radius)
                equal_radius_a_idx = find(radius_a_range == radius);
                
                deviation_results(trial_idx, radius_idx, radius_a_idx) = deviation_results(trial_idx, radius_idx, equal_radius_a_idx);
                segmentation_ratio_results(trial_idx, radius_idx, radius_a_idx) = segmentation_ratio_results(trial_idx, radius_idx, equal_radius_a_idx);
            else
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
                options.executeFindMotionSegments = true;
                options.executeClusterMotionSegments = false;
                options.generalized_radius_activities = radius_a;
                options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
                options.step_over_frame_penalty = options.radius;
                options.radius_clustering = options.radius;

                [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);

                frame_count = mot.nframes;
                [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
                [ deviation, segmentation_ratio, activity_count ] = tw_relative_segment_deviation(cuts, subcuts, frame_count);

                deviation_results(trial_idx, radius_idx, radius_a_idx) = deviation;
                segmentation_ratio_results(trial_idx, radius_idx, radius_a_idx) = segmentation_ratio;
                activity_count_results(trial_idx, radius_idx, radius_a_idx) = activity_count / 30;
            end
        end
    end
end

deviation_eval = squeeze(sum(deviation_results, 1)) / numel(trial_range);
segmentation_ratio_eval = squeeze(sum(segmentation_ratio_results, 1)) / numel(trial_range);
activity_count_eval = squeeze(sum(activity_count_results, 1)) / numel(trial_range);

%% visualization


contour_range = 0:0.01:0.99;

figure(1);
caxis_value = [0 0.5];
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
ax = gca;
[C,h] = contourf(radius_a_range, radius_range, deviation_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', radius_range);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Radius');
[void, min_d_index] = min(deviation_eval(:));
[min_d_row, min_d_column] = ind2sub(size(deviation_eval), min_d_index);
min_d_radius = radius_range(min_d_row);
min_d_radius_a = radius_a_range(min_d_column);
title(strcat('Deviation{ }', feature_set_name));


figure(2);
caxis_value = [0.5 1];
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
ax = gca;
[C,h] = contourf(radius_a_range, radius_range, segmentation_ratio_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', radius_range);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Radius');
[void, max_s_index] = max(segmentation_ratio_eval(:));
[max_s_row, max_s_column] = ind2sub(size(segmentation_ratio_eval), max_s_index);
max_s_radius = radius_range(max_s_row);
max_s_radius_a = radius_a_range(max_s_column);
title(strcat('Segmentation Ratio{ }', feature_set_name));

figure(3);
caxis_value = [0.5 1];
colorbar_ticks = caxis_value(1):0.02:caxis_value(2);
ax = gca;
[C,h] = contourf(radius_a_range, radius_range, activity_count_eval, contour_range);
clabel(C,h);
set(ax, 'XTick', radius_a_range);
set(ax, 'YTick', radius_range);
caxis(caxis_value);
colorbar('YTick', colorbar_ticks);
grid on;
xlabel('Radius Activities');
ylabel('Radius');
[void, max_a_index] = max(activity_count_eval(:));
[max_a_row, max_a_column] = ind2sub(size(activity_count_eval), max_a_index);
max_a_radius = radius_range(max_a_row);
max_a_radius_a = radius_a_range(max_a_column);
title(strcat('Relative Activity Count{ }', feature_set_name));


fprintf('Deviation Minimum := (Radius: %f, Radius Activities: %f, Value: %f)\n', min_d_radius_a, min_d_radius, deviation_eval(min_d_row, min_d_column));
fprintf('Segmentation Ratio Maximum := (Radius: %f, Radius Activities: %f, Value: %f)\n', max_s_radius_a, max_s_radius, segmentation_ratio_eval(max_s_row, max_s_column));
fprintf('Relative Activity Count Maximum := (Radius: %f, Radius Activities: %f, Value: %f)\n', max_a_radius_a, max_a_radius, activity_count_eval(max_a_row, max_a_column));