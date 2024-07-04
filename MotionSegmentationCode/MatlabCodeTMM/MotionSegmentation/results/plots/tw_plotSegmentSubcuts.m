function [] = tw_plotSegmentSubcuts( ax, subcuts )
%TW_PLOTSEGMENTSUBCUTS Summary of this function goes here
%   Detailed explanation goes here

    ylim = get(ax,'ylim');
    for i=1:numel(subcuts)
        line([subcuts(i) subcuts(i)], ylim, 'color', 'k', 'linewidth', 1);
    end
end

