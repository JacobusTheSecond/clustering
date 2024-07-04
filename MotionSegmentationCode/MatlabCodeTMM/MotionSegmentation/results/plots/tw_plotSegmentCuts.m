function [] = tw_plotSegmentCuts( ax, cuts )
%TW_PLOTSEGMENTCUTS Summary of this function goes here
%   Detailed explanation goes here

    ylim = get(ax,'ylim');
    for i=1:numel(cuts)
        line([cuts(i) cuts(i)], ylim, 'color', 'k', 'linewidth', 2);
    end
end
