function [ segment_start, segment_end, is_transition, gt_segment_color ] = tw_gt_to_segments( gt_data, annotation_data )
%TW_GT_TO_SEGMENTS Summary of this function goes here
%   Detailed explanation goes here
    
    gt_data_count = size(gt_data, 1);
    segment_count = gt_data_count * 2 - 1;

    segment_start = nan(segment_count, 1);
    segment_end = nan(segment_count, 1);
    is_transition = false(segment_count, 1);
    gt_segment_color = nan(segment_count, 3);
    
    transition_color = [0, 0, 0];
    unique_annotations = unique(annotation_data, 'stable');
    segment_colors = tw_segment_colors();

    start = 1;
    segment = 1;
    for i=1:gt_data_count
        segment_ends = gt_data(i, :);
        min_end = min(segment_ends);
        max_end = max(segment_ends);
        
        % store segment
        segment_start(segment) = start;
        segment_end(segment) = min_end;
        is_transition(segment) = false;
        [~, color_idx] = ismember(annotation_data{i}, unique_annotations);
        gt_segment_color(segment, :) = segment_colors(color_idx, :);
        
        if i ~= gt_data_count
            % store transition
            segment_start(segment+1) = min_end;
            segment_end(segment+1) = max_end;
            is_transition(segment+1) = true;
            gt_segment_color(segment+1, :) = transition_color;
        end
        
        start = max_end;
        segment = segment + 2;
    end
end

