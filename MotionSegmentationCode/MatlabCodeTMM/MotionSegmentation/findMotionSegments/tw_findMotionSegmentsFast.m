function [ subcuts,actlens,actts ] = tw_findMotionSegmentsFast( cuts, options, nnidx, nndists, is_mirrored )
%TW_FINDMOTIONSEGMENTSFAST Summary of this function goes here
%   Detailed explanation goes here
    
    DEBUG_FIND_MOTION_SEGMENTS = options.DEBUG_FIND_MOTION_SEGMENTS;
    allowedSteps = options.allowedSteps;
    step_over_frame_penalty = options.step_over_frame_penalty;
    min_path_length = options.min_path_length;
    symmetry_coverage = options.symmetry_coverage;
    min_slope = options.min_slope_motion_segments;
    max_slope = options.max_slope_motion_segments;
    min_cut_distance = options.min_cut_distance;
    use_fast_dag = options.use_fast_dag;
    use_create_dag_per_range = options.use_create_dag_per_range;
    use_old_motion_segmentation = options.use_old_motion_segmentation;

    frame_count = size(nnidx, 2);

    nnidx_filtered = nnidx;
    nndists_filtered = nndists;
    
    if (DEBUG_FIND_MOTION_SEGMENTS)
        if (is_mirrored)
            figure_offset = 5;
            title_ext = ' Mirrored';
        else
            figure_offset = 0;
            title_ext = '';
        end
            
        debug_subcuts = [];
        figure(31 + figure_offset); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title(strcat('Paths', title_ext));
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(32 + figure_offset); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title(strcat('DAG Regions', title_ext));
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(33 + figure_offset); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title(strcat('Paths, Subcuts and excluded Subcuts', title_ext));
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(1);
        figure(34 + figure_offset); 
        tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
        title(strcat('Paths and Subcuts', title_ext));
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
    
    subcuts = [];
    
    activity_count = numel(cuts) + 1;
    
    actlens = zeros(1,activity_count);
    actts = zeros(1,activity_count);
    
    % for all activities
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
        
        actlens(activity) = ef-sf+1;
        t1 = toc;
              
        is_cut_to_range = true;
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
        
        if (DEBUG_FIND_MOTION_SEGMENTS)
            % find all unfiltered subcuts
            debug_paths_subcuts = zeros(1, numel(paths));
            for debug_i = 1:numel(paths)
                debug_path = paths{debug_i};
                % subcut is at the reference start node of the path
                debug_subcut = debug_path(end, 2);
                debug_paths_subcuts(debug_i) = debug_subcut;
            end
            debug_subcuts = [debug_subcuts debug_paths_subcuts];
            
            % find all invalid paths
            [debug_paths, debug_dists, debug_is_valid_symmetry] = tw_checkSoftSymmetry(paths, dists, symmetry_coverage);
            debug_is_valid_path = debug_is_valid_symmetry;
            debug_is_valid = false(1, numel(debug_paths));
            
            for debug_i = 1:numel(debug_paths)
                debug_path = debug_paths{debug_i};

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
            
            debug_is_valid_path(debug_is_valid_path) = debug_is_valid;
            
            figure(31 + figure_offset);
            tw_plotPathsColor(paths(debug_is_valid_path), 'r');
            tw_plotPathsColor(paths(~debug_is_valid_path), 'y');
            figure(32 + figure_offset);
            tw_plotPathsColor(paths(debug_is_valid_path), 'r');
            tw_plotPathsColor(paths(~debug_is_valid_path), 'y');
            figure(33 + figure_offset);
            tw_plotPathsColor(paths(debug_is_valid_path), 'r');
            tw_plotPathsColor(paths(~debug_is_valid_path), 'y');
            figure(34 + figure_offset);
            tw_plotPathsColor(paths(debug_is_valid_path), 'r');
            figure(1);
        end
            
        % Check for soft symmetry: Any valid path except for first diagonal
        % should have a pseudo-symmetrical twin (use first diagonal as symmtry axis)
        [paths, dists, ~] = tw_checkSoftSymmetry(paths, dists, symmetry_coverage);
        
        subcuts_activity = zeros(1, numel(paths));
        path_lengths = zeros(1, numel(paths));
        
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
                % subcut is at the reference start node of the path
                subcut = path(end, 2);
                
                % check if subcut is too close to the cuts
                if subcut >= (sf + min_cut_distance) && subcut <= (ef - min_cut_distance)
                    path_length = size(path, 1);

                    subcuts_activity(j) = subcut;
                    path_lengths(j) = path_length;
                end
            end
        end
        
        is_valid_subcut = subcuts_activity ~= 0;
        subcuts_activity = subcuts_activity(is_valid_subcut);
        path_lengths = path_lengths(is_valid_subcut);
        
        subcuts_activity = tw_checkSubcutDistance(subcuts_activity, path_lengths, min_cut_distance);
        
        t2 = toc;
        actts(activity) = t2-t1;
        
        % sort subcuts
        subcuts_activity = sort(subcuts_activity);
        
        subcuts = [subcuts subcuts_activity];
    end
    
    if (DEBUG_FIND_MOTION_SEGMENTS)
        figure(31 + figure_offset);
        tw_plotCuts(cuts);
        figure(32 + figure_offset);
        tw_plotActivitySquares(cuts, frame_count);
        tw_plotCuts(cuts);
        figure(33 + figure_offset);
        [debug_is_valid] = ismember(debug_subcuts, subcuts);
        debug_filtered = debug_subcuts(~debug_is_valid);
        debug_valid = debug_subcuts(debug_is_valid);
        tw_plotSubCutsColor(debug_filtered, 'k');
        tw_plotSubCutsColor(debug_valid, 'w');
        tw_plotCuts(cuts);
        figure(34 + figure_offset);
        tw_plotSubCuts(subcuts);
        tw_plotCuts(cuts);
        figure(1);
    end
end

function [ valid_subcuts ] = tw_checkSubcutDistance(subcuts, path_lengths, min_cut_distance)
    valid_subcuts = zeros(1, numel(subcuts));
    
    [ path_lengths, sort_idx ] = sort(path_lengths, 'descend');
    subcuts = subcuts(sort_idx);
    % It is preferable to keep the subcuts, which have been determined by a
    % long path. A long path corresponds to a large connected component and
    % subcuts from large CCs are more important then subcuts from small CCs.

    for i = 1:numel(path_lengths)
       if path_lengths(i) == 0
           continue
       end
       
       valid_subcuts(i) = subcuts(i);
       
       % for all smaller paths
       for j = i+1:numel(path_lengths)
           if path_lengths(j) == 0
               continue
           end
           
           if subcuts(j) > subcuts(i)-min_cut_distance && subcuts(j) < subcuts(i)+min_cut_distance
               path_lengths(j) = 0;
           end
       end
    end
    
    valid_subcuts = valid_subcuts(valid_subcuts > 0);
end