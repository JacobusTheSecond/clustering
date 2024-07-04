function [ deviation, segmentation_ratio, activity_count, segment_count ] = tw_relative_segment_deviation( cuts, subcuts, frame_count )
%TW_RELATIVE_SEGMENT_DEVIATION Summary of this function goes here
%   Detailed explanation goes here

    allcuts = sort([cuts subcuts]);
    allcuts_plus_border = [1 allcuts frame_count];
    segment_count = numel(allcuts_plus_border) - 1;

    if (numel(subcuts) == 0)
        deviation = 0;
        segmentation_ratio = 0;
        return 
    end

    [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
    activity_count = numel(cut_segment_start);
    subcuts_count = numel(subcuts);
    subcuts_start_idx = 1;
    subcuts_end_idx = 1;
    
    activity_deviations = nan(1, activity_count);
    segmented_frames = zeros(1, activity_count);
    
    for activity=1:activity_count
        cut_start = cut_segment_start(activity);
        cut_end = cut_segment_end(activity);

        subcuts_start = subcuts(subcuts_start_idx);
        if ~(subcuts_start >= cut_start && subcuts_start < cut_end)
            continue;
        end
        
        subcuts_idx = subcuts_start_idx;
        while (subcuts_idx <= subcuts_count && subcuts(subcuts_idx) < cut_end)
            subcuts_end_idx = subcuts_idx;
            subcuts_idx = subcuts_idx + 1;
        end

        idx_range = subcuts_start_idx:subcuts_end_idx;
        if (numel(idx_range) >= 2)
            subcuts_activity = subcuts(idx_range);
            all_cuts = horzcat(cut_start, subcuts_activity, cut_end);
            segment_lengths = diff(all_cuts);
            median_length = median(segment_lengths);
            deviation_lengths = abs(segment_lengths - median_length);
            relative_deviations = deviation_lengths / median_length;
            activity_deviation = sum(relative_deviations) / numel(relative_deviations);
            
            activity_deviations(activity) = activity_deviation;
            segmented_frames(activity) = cut_end - cut_start;
        end

        % set subcuts_start_idx for the next iteration
        if (subcuts_end_idx < subcuts_count)
            subcuts_start_idx = subcuts_end_idx + 1;
        end
    end
    
    is_valid = ~isnan(activity_deviations);
    activity_deviations = activity_deviations(is_valid);
    deviation = sum(activity_deviations) / numel(activity_deviations);
        
    %segmentation_ratio = sum(is_valid) / activity_count;
    segmentation_ratio = (sum(segmented_frames) + 1) / frame_count;
end
