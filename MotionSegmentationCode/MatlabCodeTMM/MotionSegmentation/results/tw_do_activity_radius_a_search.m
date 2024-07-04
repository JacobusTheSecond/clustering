clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

trial_range = 1:14;

% ref_radius_a = 35;
% offset_radius_a = 24;
% step_radius_a = 3;
ref_radius_a = 30;
offset_radius_a = 12;
step_radius_a = 0.5;
radius_a_range = (ref_radius_a - offset_radius_a):step_radius_a:(ref_radius_a + offset_radius_a);
radius = max(radius_a_range);

[ gt, annotation ] = tw_gt_86();

strict_results = zeros(numel(trial_range), numel(radius_a_range));
tolerant_results = zeros(numel(trial_range), numel(radius_a_range));

for trial_idx=1:numel(trial_range)
    trial = trial_range(trial_idx);
    [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );
    
    gt_trial = gt{trial};
    annotation_trial = annotation{trial};
    [ gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color ] = tw_gt_to_segments(gt_trial, annotation_trial);
    
    pre_options.feature_set = 'e15_flex';
    pre_options.frame_offsets = [-10 -5 0 5 10];
    pre_options.use_mirror_motion = false;
    pre_options.use_feature_projection = false;
    pre_options.feature_projection_k = 32;
    pre_options.use_normal_weights = false;
    pre_options.normal_weights_sigma = 4;
    pre_options.k = 10000;%600;
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
        
        strict_results(trial_idx, radius_idx) = strict;
        tolerant_results(trial_idx, radius_idx) = tolerant;
    end
end

strict_eval = sum(strict_results, 1) / numel(trial_range);
tolerant_eval = sum(tolerant_results, 1) / numel(trial_range);

%% visualization

figure(1);
plot(radius_a_range, strict_eval, radius_a_range, tolerant_eval);
xlabel('Radius Activities');
ylabel('Accuracy');
legend('Strict Evaluation', 'Tolerant Evaluation');
ax = gca;
set(ax, 'XTick', radius_a_range);
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'on');