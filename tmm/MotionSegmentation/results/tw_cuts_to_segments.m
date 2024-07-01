function [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count )
%TW_CUTS_TO_SEGMENTS Summary of this function goes here
%   Detailed explanation goes here

    cut_segment_count = numel(cuts) + 1;

    cut_segment_start = nan(cut_segment_count, 1);
    cut_segment_end = nan(cut_segment_count, 1);

    % for all activities
    for activity = 1:cut_segment_count
        if activity == 1
            sf = 1;
        else
            sf = cuts(activity - 1);
        end
        
        if activity == cut_segment_count
            ef = frame_count;
        else
            ef = cuts(activity);
        end
        
        cut_segment_start(activity) = sf;
        cut_segment_end(activity) = ef;
    end
end
