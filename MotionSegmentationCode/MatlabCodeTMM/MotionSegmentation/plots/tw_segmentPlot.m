function [] = tw_segmentPlot( cuts, subcuts, submots, frame_count )
%TW_TIMINGPLOT Summary of this function goes here
%   Detailed explanation goes here

    allcuts = sort([cuts subcuts]);
    segment_count = numel(allcuts) + 1;
    activity_count = numel(cuts) + 1;
    permColorIndexes = randperm(activity_count);
    colors = jet(activity_count);
    colors = colors(permColorIndexes, :);
    color_idx = 1;
    
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
        
        w  = ef-sf;
        if ismember(sf, cuts)
            color_idx = color_idx + 1;
        end
        color = colors(color_idx, :);
           
        rectangle('Position', [sf, 0 , w, 1], 'Facecolor', color, 'linewidth', 2);
    end

    set(gca,'xlim',[1 frame_count]);
    %set(gcf,'position',[100 100 1800 450])

    %cellfun(@tw_timingPlotRectangle, submots);
end

% function [submot] = tw_timingPlotRectangle( submot )
%     sf = submot.sf;
%     ef = submot.ef;
%     color = submot.color;
%     
%     w  = ef-sf;
%     rectangle('Position', [sf, 0 , w, 1], 'Facecolor', color, 'linewidth', 2);
% end