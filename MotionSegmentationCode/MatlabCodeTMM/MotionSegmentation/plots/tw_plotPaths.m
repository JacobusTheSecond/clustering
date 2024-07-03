function [] = tw_plotPaths(paths)    
    for i=1:numel(paths)
       path = paths{i};
       line(path(:, 2), path(:, 1),'color','red','linewidth',1);
    end
end
