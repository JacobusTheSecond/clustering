function [ mots, submots, comps ] = tw_clusterMotionSegments( cuts, subcuts, mot, options, nnidx, nndists )
%TW_CLUSTERMOTIONSEGMENTS Summary of this function goes here
%   Detailed explanation goes here

    filter_radius = options.radius_clustering;
    allowedSteps = options.allowedSteps;
    step_over_frame_penalty = options.step_over_frame_penalty;
    min_path_length = options.min_path_length;
    min_slope = options.min_slope_clustering;
    max_slope = options.max_slope_clustering;
    coverage_clustering = options.coverage_clustering;
    use_fast_dag = options.use_fast_dag;
    use_create_dag_per_range = options.use_create_dag_per_range;
    use_old_motion_segmentation = options.use_old_motion_segmentation;
    DEBUG_CLUSTER_MOTION_SEGMENTS = options.DEBUG_CLUSTER_MOTION_SEGMENTS;
    
    frame_count = size(nnidx, 2);

    allcuts = sort([cuts subcuts]);
    allcuts_plus_border = [1 allcuts frame_count];

    segment_count = numel(allcuts) + 1;
    segment_adjacency_matrix  = zeros(segment_count);
    
    [nnidx_filtered nndists_filtered] = tw_filterRadius(filter_radius, nnidx, nndists);
    
    if (DEBUG_CLUSTER_MOTION_SEGMENTS)
        figure(41); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title('Paths');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(42); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title('Paths and Clusters');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(1);
    end  
    
    if (~use_create_dag_per_range)
        if (use_fast_dag)
            A = tw_buildDagFast(nnidx_filtered, nndists_filtered, allowedSteps, step_over_frame_penalty, frame_count);
        else
            A = tw_buildDAGmatrix(nnidx_filtered, nndists_filtered, allowedSteps, step_over_frame_penalty);
        end
    end
        
    % for all segments
    for segment = 1:segment_count
        if segment == 1
            sf = 1;
        else
            sf = allcuts(segment - 1);
        end

        if segment == segment_count
            ef = frame_count;
        else
            ef = allcuts(segment);
        end
        
        %consider strip of Afull corresponding to current segment and look for
        %suitable paths

        is_cut_to_range = false;
        if (use_create_dag_per_range)
            [ A_cut, nnidx_cut, nndists_cut ] = tw_createCutDagSingleStartEnd( sf, ef, is_cut_to_range, nnidx_filtered, nndists_filtered, allowedSteps, step_over_frame_penalty );
        else
            [ A_cut, nnidx_cut, nndists_cut ] = tw_cutDagSingleStartEnd( A, sf, ef, is_cut_to_range, nnidx_filtered, nndists_filtered );
        end
        
        if (use_old_motion_segmentation)
            [ paths, dists ] = tw_getSubsegments(nnidx_cut, nndists_cut, A_cut, sf);
        else
            [ paths, dists ] = tw_fastSubsegments(A_cut, sf, min_path_length, nnidx_cut, nndists_cut);
        end
        
        if (DEBUG_CLUSTER_MOTION_SEGMENTS)
            % find all invalid paths
            debug_is_valid = false(1, numel(paths));
            
            for debug_i = 1:numel(paths)
                debug_path = paths{debug_i};

                % a path ends at index 1 and starts at index end
                debug_path_range = debug_path(1, 1) - debug_path(end, 1);
                debug_ref_path_range = debug_path(1, 2) - debug_path(end, 2);

                if debug_path_range == 0
                    continue
                end

                debug_slope = debug_ref_path_range / debug_path_range;

                if debug_slope >= min_slope && debug_slope <= max_slope
                    debug_is_valid(debug_i) = true;
                end
            end
                        
            figure(41);
            tw_plotPathsColor(paths(debug_is_valid), 'r');
            tw_plotPathsColor(paths(~debug_is_valid), 'y');
            figure(42);
            tw_plotPathsColor(paths(debug_is_valid), 'r');
            tw_plotPathsColor(paths(~debug_is_valid), 'y');
            figure(1);
        end
            
        segment_frames = sf:ef;

        % for all segment paths
        for j = 1:numel(paths)
            path = paths{j};
                        
            % a path ends at index 1 and starts at index end
            path_range = path(1, 1) - path(end, 1);
            ref_path_range = path(1, 2) - path(end, 2);
            
            if path_range == 0
                continue
            end
            
            slope = ref_path_range / path_range;
            
            if slope >= min_slope && slope <= max_slope
                path_length = size(path, 1);
                %path_frames = path(end, 1):path(1, 1);
                ref_path_frames = path(end, 2):path(1, 2);
                
                first_ref_frame = ref_path_frames(1);
                last_ref_frame = ref_path_frames(end);
                
                %first_ref_closest_subcut_idx = tw_getClosestFrameIndex(first_ref_frame, allcuts_plus_border);
                %last_ref_closest_subcut_idx = tw_getClosestFrameIndex(last_ref_frame, allcuts_plus_border);
                
                first_ref_frame_closest_allcut_idx = tw_getClosestFrameListIndexLeft(first_ref_frame, allcuts_plus_border);
                last_ref_frame_closest_allcut_idx = tw_getClosestFrameListIndexRight(last_ref_frame, allcuts_plus_border);
                
                allcuts_idxs = first_ref_frame_closest_allcut_idx:last_ref_frame_closest_allcut_idx;
                            
                path_flipped = flipud(path);

                
                %TODO remove
%                 if (j == 8 && segment == 2)
%                     dbstop in title
%                     title('abc')
%                 end
%                 if (j == 4 && segment == 2)
%                     dbstop in title
%                     title('abc')
%                 end
                
                
                % iterate over all ref segments
                for r = 1:numel(allcuts_idxs)-1
                    % get frames between two subcuts (segment frames)
                    ref_segment_frames = allcuts_plus_border( allcuts_idxs(r) ):allcuts_plus_border( allcuts_idxs(r+1) );
                    [path_start_left, path_start_right] = tw_getClosestListIndices(ref_segment_frames(1), path_flipped(:, 2), path_length);
                    [path_end_left, path_end_right] = tw_getClosestListIndices(ref_segment_frames(end), path_flipped(:, 2), path_length);
                    path_start = round( ( path_flipped(path_start_left, 1) + path_flipped(path_start_right, 1) ) / 2 );
                    path_end = round( ( path_flipped(path_end_left, 1) + path_flipped(path_end_right, 1) ) / 2 );
                    path_frames = path_start:path_end;
                    
                    if getcoverage_local(ref_path_frames, path_frames, ref_segment_frames, segment_frames, coverage_clustering)
                        % Since r <= (allcuts last index - 1), allcuts_idxs(r) 
                        % directly corresponds to the referenced segment number
                        ref_segment = allcuts_idxs(r);              
                        segment_adjacency_matrix(segment, ref_segment) = 1;
                    end
                end
            end
        end
    end

    % clmat logical matrix
    % The rows represent the subcut segments and the columns the segments
    % that are referenced by the subcut's paths.
    % (referenced = second path parameter nnidx of frame)
    % The last row is always filled with zeros, since the size of clmat is
    % [numel(allsubcuts) + 1] X [numel(allsubcuts) + 1].
    % The last column however can contain ones.
    segment_adjacency_matrix = sparse(segment_adjacency_matrix);
    % Describes to which CC each segment belongs
    comps = components_mex(segment_adjacency_matrix);

    activity_count = numel(cuts) + 1;
    mots = cell(activity_count, 1);
    
    for activity = 1:activity_count
        if activity == 1
            sf = 1;
        else
            sf = cuts(activity - 1);
        end

        if activity == activity_count
            ef = frame_count;
        else
            ef = cuts(activity);
        end

        mots{activity} = cutMotion(mot, sf, ef);
    end
    
    comps_count = max(comps);
    colors = jet(comps_count);
    submots = cell(segment_count, 1);
    
    for segment = 1:segment_count
        if segment == 1
            sf = 1;
        else
            sf = allcuts(segment - 1);
        end

        if segment == segment_count
            ef = frame_count;
        else
            ef = allcuts(segment);
        end

        submots{segment} = cutMotion(mot, sf, ef);
        % Get RGB color for submot segment. The color indicates the
        % connected component the segment belongs to.
        submots{segment}.color = colors(comps(segment), :);
        submots{segment}.sf = sf;
        submots{segment}.ef = ef;
    end
    
    if (DEBUG_CLUSTER_MOTION_SEGMENTS)
        figure(41);
        tw_plotSubCuts(subcuts);
        tw_plotCuts(cuts);
        figure(42);
        tw_plotClusterSquares(cuts, subcuts, comps, frame_count);
        tw_plotSubCuts(subcuts);
        tw_plotCuts(cuts);
        figure(1);
    end
end

function coverage = getcoverage_local(ref_path_frames, path_frames, ref_segment_frames, segment_frames, coverage_clustering)
    coverage = 1;
    
    xcov = sum( ismember(ref_segment_frames, ref_path_frames) ) / numel(ref_segment_frames);
    ycov = sum( ismember(segment_frames, path_frames) ) / numel(segment_frames);

    if xcov < coverage_clustering
        coverage = 0;
    elseif ycov < coverage_clustering
        coverage = 0;
    end
end

function [closest_index_left, closest_index_right] = tw_getClosestListIndices(item, list, max_index)
    distance_list = list - item;
    closest_index_left = sum(distance_list <= 0);
    closest_index_left = max(closest_index_left, 1);
    if (distance_list(closest_index_left) < 0)
        closest_index_right = min(closest_index_left + 1, max_index);
    else
        closest_index_right = closest_index_left;
    end
end

function closest_index_left = tw_getClosestFrameListIndexLeft(frame, frame_list)
    distance_list = frame_list - frame;
    % sum(distance_list <= 0) ensures that if the given frame is equal 
    % to a frame in the frame list, then the index of the frame of 
    % the frame list is returned. Otherwise the index of the
    % closest frame in frame list at the left side is returned.
    closest_index_left = sum(distance_list <= 0);
end

function closest_index_right = tw_getClosestFrameListIndexRight(frame, frame_list)
    distance_list = frame_list - frame;
    % sum(distance_list < 0) + 1 ensures that if the given frame is equal 
    % to a frame in the frame list, then the index of the frame of 
    % the frame list is returned. Otherwise the index of the
    % closest frame in frame list at the right side is returned.
    closest_index_right = sum(distance_list < 0) + 1;
end