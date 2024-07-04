function [ strict, tolerant, cut_segment_to_gt_segment, cut_segment_color ] = tw_segment_match( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, cut_segment_start, cut_segment_end)
%TW_SEGMENT_MATCH Summary of this function goes here
%   Detailed explanation goes here

    gt_segment_count = numel(gt_segment_start);
    cut_segment_count = numel(cut_segment_start);
    
    cut_segment_to_gt_segment = nan(cut_segment_count, 1);    
    cut_segment_color = nan(cut_segment_count, 3);
    
    %% assign cut segments to gt segments
    for cut_segment=1:cut_segment_count
        coverage = 0;
        matching_gt_segment = nan;
        
        for gt_segment=1:gt_segment_count
            gt_start = gt_segment_start(gt_segment);
            gt_end = gt_segment_end(gt_segment);
            cut_start = cut_segment_start(cut_segment);
            cut_end = cut_segment_end(cut_segment);

            intersect_start = max(gt_start, cut_start);
            intersect_end = min(gt_end, cut_end);

            intersect_range = max(intersect_end - intersect_start + 1, 0);
            cut_range = cut_end - cut_start + 1;
            
            coverage_percent = intersect_range / cut_range;
            
            if coverage_percent > coverage
                matching_gt_segment = gt_segment;
                coverage = coverage_percent;
            end
        end
        
        cut_segment_to_gt_segment(cut_segment) = matching_gt_segment;
        cut_segment_color(cut_segment, :) = gt_segment_color(matching_gt_segment, :);
    end
    
    %% compute strict and tolerant value
    strict_intersect = 0;
    strict_count = 0;
    tolerant_intersect = 0;
    tolerant_count = 0;
    
    for gt_segment=1:gt_segment_count
        gt_start = gt_segment_start(gt_segment);
        gt_end = gt_segment_end(gt_segment);
        gt_range = gt_end - gt_start + 1;
        
        gt_cuts = cut_segment_to_gt_segment == gt_segment;
        
        cuts_start = min(cut_segment_start(gt_cuts));
        cuts_end = max(cut_segment_end(gt_cuts));
        %cuts_range = cuts_end - cuts_start + 1;
        
        if isempty(cuts_start)
            intersect_range = 0;
        else
            intersect_start = max(gt_start, cuts_start);
            intersect_end = min(gt_end, cuts_end);
            intersect_range = intersect_end - intersect_start + 1;
        end
        
        strict_intersect = strict_intersect + intersect_range;
        strict_count = strict_count + gt_range;
            
        if ~gt_is_transition(gt_segment)
            tolerant_intersect = tolerant_intersect + intersect_range;
            tolerant_count = tolerant_count + gt_range;
        end
    end
    
    strict = strict_intersect / strict_count;
    tolerant = tolerant_intersect / tolerant_count;   
end

