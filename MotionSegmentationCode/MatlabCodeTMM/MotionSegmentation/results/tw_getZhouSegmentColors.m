function [ zhou_segment_color ] = tw_getZhouSegmentColors( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, zhou_segment_start, zhou_segment_end, zhou_clustering )
%TW_GETZHOUSEGMENTCOLORS Summary of this function goes here
%   Detailed explanation goes here

    gt_segment_count = numel(gt_segment_start);
    zhou_segment_count = numel(zhou_segment_start);
    
    zhou_segment_color = zeros(zhou_segment_count, 3);
    
    zhou_cluster_start = nan(1, zhou_segment_count+1);
    zhou_cluster_end = nan(1, zhou_segment_count+1);
    zhou_cluster_id = nan(1, zhou_segment_count+1);
    zhou_cluster_start_idx = nan(1, zhou_segment_count+1);
    zhou_cluster_end_idx = nan(1, zhou_segment_count+1);
    
    start_idx = 1;
    end_idx = 1;
    zhou_cluster = zhou_clustering(1);
    for i=1:zhou_segment_count
        if (zhou_clustering(i) == zhou_cluster)
            end_idx = i;
        else
            zhou_cluster_start(i) = zhou_segment_start(start_idx);
            zhou_cluster_end(i) = zhou_segment_end(end_idx);
            zhou_cluster_id(i) = zhou_cluster;
            zhou_cluster_start_idx(i) = start_idx;
            zhou_cluster_end_idx(i) = end_idx;
            
            start_idx = i;
            end_idx = i;
            zhou_cluster = zhou_clustering(i);    
        end
    end
    
    zhou_cluster_start(end) = zhou_segment_start(start_idx);
    zhou_cluster_end(end) = zhou_segment_end(end_idx);
    zhou_cluster_id(end) = zhou_cluster;
    zhou_cluster_start_idx(end) = start_idx;
    zhou_cluster_end_idx(end) = end_idx;
    
    is_valid = ~isnan(zhou_cluster_start);
    zhou_cluster_start = zhou_cluster_start(is_valid);
    zhou_cluster_end = zhou_cluster_end(is_valid);
    zhou_cluster_id = zhou_cluster_id(is_valid);
    zhou_cluster_start_idx = zhou_cluster_start_idx(is_valid);
    zhou_cluster_end_idx = zhou_cluster_end_idx(is_valid);
    
    zhou_cluster_count = numel(zhou_cluster_start);
    gt_segment_to_cluster_id = nan(1, gt_segment_count);
    gt_segment_coverage = zeros(1, gt_segment_count);
    cluster_id_to_gt_segment = nan(1, zhou_cluster_count);
    
    for zhou_cluster = 1:zhou_cluster_count
        cluster_coverage = 0;
        matching_gt_segment = nan;
        
        cluster_start = zhou_cluster_start(zhou_cluster);
        cluster_end = zhou_cluster_end(zhou_cluster);
        cluster_id = zhou_cluster_id(zhou_cluster);
        cluster_start_idx = zhou_cluster_start_idx(zhou_cluster);
        cluster_end_idx = zhou_cluster_end_idx(zhou_cluster);
        
        if ( isnan( cluster_id_to_gt_segment(cluster_id) ) )     
            for gt_segment = 1:gt_segment_count
                gt_start = gt_segment_start(gt_segment);
                gt_end = gt_segment_end(gt_segment);
                gt_transition = gt_is_transition(gt_segment);

                if (~gt_transition)
                    intersect_start = max(gt_start, cluster_start);
                    intersect_end = min(gt_end, cluster_end);

                    intersect_range = max(intersect_end - intersect_start + 1, 0);
                    cluster_range = cluster_end - cluster_start + 1;

                    coverage_percent = intersect_range / cluster_range;

                    if coverage_percent > cluster_coverage
                        matching_gt_segment = gt_segment;
                        cluster_coverage = coverage_percent;
                        gt_coverage = intersect_range / (gt_end - gt_start + 1);
                    end
                end
            end

            if isnan(gt_segment_to_cluster_id(matching_gt_segment))
                gt_segment_to_cluster_id(matching_gt_segment) = cluster_id;
                cluster_id_to_gt_segment(cluster_id) = matching_gt_segment;
                gt_segment_coverage(matching_gt_segment) = gt_coverage;
            elseif gt_segment_coverage(matching_gt_segment) < gt_coverage 
                old_cluster_id = gt_segment_to_cluster_id(matching_gt_segment);
                cluster_id_to_gt_segment(old_cluster_id) = nan;
                
                gt_segment_to_cluster_id(matching_gt_segment) = cluster_id;
                cluster_id_to_gt_segment(cluster_id) = matching_gt_segment;
                gt_segment_coverage(matching_gt_segment) = gt_coverage;
            end
        end
    end
    
    valid_cluster_ids = find(~isnan(cluster_id_to_gt_segment));
    invalid_cluster_ids = find(isnan(cluster_id_to_gt_segment));
    
    for i=1:numel(valid_cluster_ids)
       cluster_id = valid_cluster_ids(i);
       cluster_idx = (zhou_clustering == cluster_id);
       segment_color = gt_segment_color(cluster_id_to_gt_segment(cluster_id), :);
       zhou_segment_color(cluster_idx, :) = repmat( segment_color, sum(cluster_idx), 1);
    end
    
    for i=1:numel(invalid_cluster_ids)
        cluster_id = invalid_cluster_ids(i);
        cluster_idx = (zhou_clustering == cluster_id);
        segment_color = [0 0 0];
        zhou_segment_color(cluster_idx, :) = repmat( segment_color, sum(cluster_idx), 1);
    end
end

