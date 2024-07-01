function [] = tw_plotCuts(cuts)
%TW_PLOTCUTS Summary of this function goes here
%   Detailed explanation goes here

    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    for i=1:numel(cuts)
        line([cuts(i) cuts(i)], ylim,'color','black','linewidth',2);
        line(xlim, [cuts(i) cuts(i)],'color','black','linewidth',2);
    end
end