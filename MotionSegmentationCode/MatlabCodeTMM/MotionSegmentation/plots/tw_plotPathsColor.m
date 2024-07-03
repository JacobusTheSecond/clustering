function [] = tw_plotPathsColor(paths, color)    
    for i=1:numel(paths)
       path = paths{i};
       line(path(:, 2), path(:, 1),'color',color,'linewidth',3);
    end
end
