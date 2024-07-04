function [] = tw_plotSubCuts(allsubcuts)
%TW_PLOTSUBCUTS Summary of this function goes here
%   Detailed explanation goes here

    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    for i=1:numel(allsubcuts)
        line([allsubcuts(i) allsubcuts(i)], ylim,'color',[0.5 0.5 0.5]);
        line(xlim, [allsubcuts(i) allsubcuts(i)],'color',[0.5 0.5 0.5]);
    end
end