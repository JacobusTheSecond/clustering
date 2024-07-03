function [] = tw_plotSubCutsColor(allsubcuts, color)
%TW_PLOTSUBCUTS Summary of this function goes here
%   Detailed explanation goes here

    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    for i=1:numel(allsubcuts)
        line([allsubcuts(i) allsubcuts(i)], ylim,'color',color);
        line(xlim, [allsubcuts(i) allsubcuts(i)],'color',color);
    end
end