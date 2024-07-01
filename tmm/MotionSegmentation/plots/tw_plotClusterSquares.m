function [] = tw_plotClusterSquares( cuts, subcuts, comps, frame_count )
%TW_PLOTCLUSTERSQUARES Summary of this function goes here
%   Detailed explanation goes here

    allcuts = sort([cuts subcuts]);
    allcuts_plus_border = [1 allcuts frame_count];
    comps_count = max(comps);
    color_step = 0.66 / comps_count;
    colors = 0.33+color_step:color_step:1;
    colors = fliplr(colors);
    
    for row=1:numel(allcuts_plus_border)-1
        for column=1:numel(allcuts_plus_border)-1
            if (comps(row) == comps(column))
                color = colors(comps(row));
%                 x = allcuts_plus_border(column);
%                 w = allcuts_plus_border(column+1) - x;
%                 y = allcuts_plus_border(row);
%                 h = allcuts_plus_border(row+1) - y;
%                 rectangle('Position', [x, y, w, h], 'Facecolor', color, 'Edgecolor', 'none');
                x1 = allcuts_plus_border(column);
                x2 = allcuts_plus_border(column+1);
                y1 = allcuts_plus_border(row);
                y2 = allcuts_plus_border(row+1);
                p = patch([x1 x1 x2 x2], [y1 y2 y2 y1], color);
                set(p, 'FaceAlpha', 0.75);    
            end
        end
    end
end

