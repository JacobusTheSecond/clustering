function [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings, meta] = tw_segmentation(skel, mot, varargin)
%TW_SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

    %% Global options
    options.DEBUG = true;
    options.DEBUG_FIND_MOTION_ACTIVITIES = false;
    options.DEBUG_FIND_MOTION_SEGMENTS = false;
    options.DEBUG_CLUSTER_MOTION_SEGMENTS = false;
    
    options.executeChangeFrameRate = true;
    options.executeCreateFeatureSet = true;
    options.fmat = [];
    options.fmat_mirror = [];
    options.executeFindMotionSimilarities = true;
    options.nnidx = [];
    options.nndists = [];
    options.nnidx_mirror = [];
    options.nndists_mirror = [];
    options.executeFindMotionSegments = true;
    options.executeClusterMotionSegments = true;
    
    options.frameRate = 30;
    %options.frameRate = 60;
    %options.frameRate = 120;
    
    %options.feature_set     = 'e15_45';
    options.feature_set     = 'e15_flex';
    %options.feature_set     = 'e30_flex';
    
    options.frame_offsets = [-10 -5 0 5 10];
    
    options.feature_projection_k = 32;
    
    % adjust frame offsets to sampling rate
    options.frame_offsets = options.frame_offsets * options.frameRate / 30;
    
    options.use_normal_weights = false;
    options.normal_weights_sigma = 4;
    
    options.k = 600;
    
    %options.radius   = 32;%64;%80;%42;%64;
    %options.radius_activities = 32;%50;%40;% * 0.66;
    
    % radius 48 good for projection [-10 -5 0 5 10];
    % radius 60 good without projection [-10 -5 0 5 10]
    % radius 125 good without projection [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10]
    % radius 68 good without projection [-10 -5 0 5 10] filtering step1 0.9 
    
    %options.radius = 68;%48;%60;%80;%42;%64;
    
    % generalized radius
    options.generalized_radius = 19;%16;%30.2;%40;%30.5;
    
    % calculate frame window specific radius 
    options.allowedSteps      = [ 1 1; ...
                                1 0; ...
                                0 1; ...
%                                 1 2; ...
%                                 2 1; ...
%                                 2 2; ...
%                                 5 5; ...
%                                 10 10; ...
                                ];
    
    
    options.min_path_length = options.frameRate / 2;
    
    options.use_mirror_motion = false;
    options.use_fast_dag = false;
    options.use_create_dag_per_range = true;
    
    options.use_old_activity_search = false;
    options.use_old_motion_segmentation = false;
    
    %% use a PCA projection for the feature set
    options.use_feature_projection = true;
    
    %% Find motion activities options
    options.generalized_radius_activities = 27;%36;%options.generalized_radius * 0.88;
   
    % calculate frame window specific activity radius 
    options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
    
    % minimum number of nearest neighbors in a blob
    options.nn_min_blob_size = options.frameRate * 4;%128;
    options.min_zero_deltas_cut = floor(options.frameRate / 6);%floor(mot.samplingRate / 3);
    options.filter_diagonal_width = floor(options.frameRate);%floor(mot.samplingRate / 2);%15;
    
    %% Find motion segments options
    options.symmetry_coverage = 0.5;
    options.min_cut_distance = options.frameRate / 4;
    options.min_slope_motion_segments = 0.5;
    options.max_slope_motion_segments = 1 / options.min_slope_motion_segments;
    
    %% Cluster motion segments options
    
    options.coverage_clustering = 0.66;
    options.min_slope_clustering = 0.5;
    options.max_slope_clustering = 1 / options.min_slope_clustering;
    
    %% Merge options
    if (nargin == 3)
        options = mergeOptions(varargin{1}, options);
    end

    options.radius = options.generalized_radius * sqrt( numel(options.frame_offsets) );
    options.step_over_frame_penalty = options.radius;
    options.radius_clustering = options.radius;
    
    %% Start algorithm
    
    timings = {};
    
    %% Set frame rate
    if (options.executeChangeFrameRate)
        mot = changeFrameRate(skel, mot, options.frameRate);
        options.k = min(mot.nframes, options.k);
    end
    
    %% Create Feature Set
    tic;
    if (options.executeCreateFeatureSet)
        [ fmat, fmat_mirror ] = tw_createFeatureSet( skel, mot, options );
    else
        fmat = options.fmat;
        fmat_mirror = options.fmat_mirror;
    end
    timings.Create_Feature_Set = toc;
    fprintf('Create Feature Set\t\t\tcompleted in % 2f seconds.\n', timings.Create_Feature_Set);
    
    %% Find Motion Similarities
    tic;
    if (options.executeFindMotionSimilarities)
        [ nnidx, nndists, nnidx_mirror, nndists_mirror ] = tw_findMotionSimilarities( fmat, fmat_mirror, options );
        
        if (options.DEBUG)
            figure(1);

            frame_count = mot.nframes;

            subplot(1,2,1);
            if (options.use_mirror_motion)
    %             tw_mergePlot(frame_count, options.radius, nnidx, nndists, nnidx_mirror, nndists_mirror, nnidx_merge, nndists_merge); 
                tw_mirrorPlot(frame_count, options.radius, nnidx, nndists, nnidx_mirror, nndists_mirror);
            else
                tw_selfSimilarityPlot(frame_count, options.radius, nnidx, nndists);
            end

            subplot(1,2,2);
            if (options.use_mirror_motion)
    %             tw_mergePlot(frame_count, options.radius, nnidx, nndists, nnidx_mirror, nndists_mirror, nnidx_merge, nndists_merge); 
                tw_mirrorPlot(frame_count, options.radius, nnidx, nndists, nnidx_mirror, nndists_mirror);
            else
                tw_selfSimilarityPlot(frame_count, options.radius, nnidx, nndists);
            end
        end
    else
        nnidx = options.nnidx;
        nndists = options.nndists;
        nnidx_mirror = options.nnidx_mirror;
        nndists_mirror = options.nndists_mirror;
    end
    timings.Find_Motion_Similarities = toc;
    fprintf('Find Motion Similarities\tcompleted in % 2f seconds.\n', timings.Find_Motion_Similarities);
    
    frame_count = size(nnidx, 2);
    last_row_invalid_count = sum(isnan(nnidx(end, :)));
    if last_row_invalid_count ~= frame_count
       truncated_count = frame_count - last_row_invalid_count;
       fprintf('Results %d from %d results are truncated. k=%d might be too small.\n', truncated_count, frame_count, options.k)
    end
 
    %% Find Motion Activities
    tic;
    [ cuts ] = tw_findActivities(nnidx, nndists, options);
    
    if (options.DEBUG)
        subplot(1,2,1);
        tw_plotCuts(cuts);
        
        subplot(1,2,2);
        tw_plotCuts(cuts);
    end
    timings.Find_Motion_Activities = toc;
    fprintf('Find Motion Activities\t\tcompleted in % 2f seconds.\n', timings.Find_Motion_Activities);


    %% Find Motion Segments
    tic;
    if (options.executeFindMotionSegments)
        [subcuts_main,meta.act_lengths] = tw_findMotionSegmentsFast(cuts, options, nnidx, nndists, false);
        if (options.use_mirror_motion)
            subcuts_mirror = tw_findMotionSegmentsFast(cuts, options, nnidx_mirror, nndists_mirror, true);
            subcuts = tw_mergeMirrorCuts(subcuts_main, subcuts_mirror, options.min_cut_distance);
        else
            subcuts = subcuts_main;
            subcuts_mirror = [];
        end
        
        if (options.DEBUG)
            subplot(1,2,2);
            tw_plotSubCuts(subcuts);
        end
    else
        subcuts = [];
        subcuts_main = [];
        subcuts_mirror = [];
    end
    timings.Find_Motion_Segments = toc;
    fprintf('Find Motion Segments\t\tcompleted in % 2f seconds.\n', timings.Find_Motion_Segments);
    
    %% Cluster Motion Segments
    tic;
    if (options.executeFindMotionSegments && options.executeClusterMotionSegments)
        [ mots, submots, comps ] = tw_clusterMotionSegments(cuts, subcuts, mot, options, nnidx, nndists);
    else
        mots={};
        submots={};
        comps={};
    end
    timings.Cluster_Motion_Segments = toc;
    fprintf('Cluster Motion Segments\t\tcompleted in % 2f seconds.\n', timings.Cluster_Motion_Segments);
end