clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

radius = 29.6;

use_feature_projection = false;

feature_set = 'e15_flex';
%feature_set = 'e30_flex';

% [0] offset gives highest values
frame_offsets = 0;

k = 1000000;
k_min = 0;

trial_range = 1:14;

for trial_idx=1:numel(trial_range)
    trial = trial_range(trial_idx);
    [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );

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

    nnidx_zero = nnidx;
    nnidx_zero(isnan(nnidx)) = 0;
    sum_zero = sum(nnidx_zero, 2);
    local_k_min = find(sum_zero, 1, 'last');
    
    if (local_k_min > k_min)
       k_min = local_k_min; 
    end
end
    
k_min = k_min + 1;
fprintf('k_min = %d\n', k_min);
