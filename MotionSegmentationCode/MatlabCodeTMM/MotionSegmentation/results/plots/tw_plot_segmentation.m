function [] = tw_plot_segmentation( segments_start, segments_end, segments_color )
%TW_PLOT_SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

    segments_count = numel(segments_start);
    xlim([segments_start(1) segments_end(end)]);

    % for all segments
    for segment = 1:segments_count
        segment_start = segments_start(segment);
        segment_end = segments_end(segment);
        segment_color = segments_color(segment, :);
            
        segment_width  = segment_end - segment_start;
           
        rectangle('Position', [segment_start, 0 , segment_width, 1], 'Facecolor', segment_color, 'linewidth', 1);
    end

end

