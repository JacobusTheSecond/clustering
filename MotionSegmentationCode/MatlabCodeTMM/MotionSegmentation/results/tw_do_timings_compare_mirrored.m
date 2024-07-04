clear all;
close all;

base_path = 'D:\mocap\cmu\all_asfamc\subjects\86\';
skel_file = '86.asf';
file_prefix = '86_';
file_suffix = '.amc';

trial_range = 1:14;
%trial_range = 1:3;

radius_a = 25;
radius = 29.6;
k = 698;

first_timings_results = cell(1, numel(trial_range));
first_names = cell(1, numel(trial_range));
second_timings_results = cell(1, numel(trial_range));
second_names = cell(1, numel(trial_range));

frame_counts = zeros(1, numel(trial_range));

for type=1:2
    for trial_idx=1:numel(trial_range)
        trial = trial_range(trial_idx);
        [ skel, mot, name ] = tw_get_mot_by_index( trial, base_path, skel_file, file_prefix, file_suffix );

        options = {};
        options.feature_set = 'e15_flex';
        options.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
        if type == 1
            options.use_mirror_motion = false;
        else
            options.use_mirror_motion = true;
        end
        options.use_feature_projection = false;  
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
        options.executeClusterMotionSegments = true;
        options.generalized_radius_activities = radius_a;
        options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
        options.step_over_frame_penalty = options.radius;
        options.radius_clustering = options.radius;
        
        options.use_create_dag_per_range = true;
    
        [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);
        if type == 1
            first_timings_results{trial_idx} = timings;
            frame_counts(trial_idx) = mot.nframes;
            first_names{trial_idx} = name;
        else
            second_timings_results{trial_idx} = timings;
            frame_counts(trial_idx) = mot.nframes;
            second_names{trial_idx} = strcat({'MM '}, name);
        end 
    end
end


%% visualization

figure(1);
timings_table = tw_timingsComparePlot(first_timings_results, second_timings_results, frame_counts, first_names, second_names, 'Average', 'MM Average', 'Normal Version vs. Mirrored Motions (MM) Version');