function [] = tw_plotCompareSegmentation(segments_cuts_list, segments_subcuts_list, segments_start_list, segments_end_list, segments_color_list, compare_names, title_name )
%TW_PLOTCOMPARESEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

    segment_data_count = numel(segments_start_list);

    if (segment_data_count < 2)
        error('Invalid Argument!\n');
    end
    
    offset_x = 0.015;
    
    subplot(segment_data_count, 1, 1);
    
    % set axes
    pos_ax1 = get(gca, 'Position');
    pos_ax2 = pos_ax1;
    pos_ax2(1) = pos_ax2(1) + offset_x;
    pos_ax2(2) = pos_ax2(2);
    pos_ax2(3) = pos_ax2(3) - offset_x;
    pos_ax2(4) = pos_ax2(4);
    ax1 = gca;
    ax2 = axes('Position', pos_ax2);
    axes(ax2);
    
    first_segments_start = segments_start_list{1};
    first_segments_end = segments_end_list{1};
    first_segments_color = segments_color_list{1};
    tw_plot_segmentation(first_segments_start, first_segments_end, first_segments_color);
    set(ax2, 'YTickLabel', []);
    title(title_name, 'Interpreter', 'none', 'FontSize', 12);
    
    axes(ax1);
    text(offset_x, 0.5, compare_names{1}, 'Interpreter', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 11);
    set(ax1, 'Visible', 'off');
        
    for i=2:segment_data_count
        subplot(segment_data_count, 1, i);
        
         % set axes
        pos_ax1 = get(gca, 'Position');
        pos_ax2 = pos_ax1;
        pos_ax2(1) = pos_ax2(1) + offset_x;
        pos_ax2(2) = pos_ax2(2);
        pos_ax2(3) = pos_ax2(3) - offset_x;
        pos_ax2(4) = pos_ax2(4);
        ax1 = gca;
        ax2 = axes('Position', pos_ax2);
        axes(ax2);
        
        segments_start = segments_start_list{i};
        segments_end = segments_end_list{i};
        segments_color = segments_color_list{i};
        tw_plot_segmentation(segments_start, segments_end, segments_color);
        cuts = segments_cuts_list{i};
        subcuts = segments_subcuts_list{i};
        tw_plotSegmentCuts( ax2, cuts );
        tw_plotSegmentSubcuts( ax2, subcuts );
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        
        axes(ax1);
        text(offset_x, 0.5, compare_names{i}, 'Interpreter', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 11);
        set(ax1, 'Visible', 'off');
    end
end

