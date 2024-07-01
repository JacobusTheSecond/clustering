function [] = tw_plotActivitySquares( cuts, frame_count )
%TW_PLOTACTIVITYSQUARES Summary of this function goes here
%   Detailed explanation goes here

    cuts_plus_border = [1 cuts frame_count];
    color_step = 0.66 / (numel(cuts_plus_border) - 1);
    colors = 0.33+color_step:color_step:1;
    
    for row=1:(numel(cuts_plus_border) - 1)
        column = row;
        color = colors(row);
                
        x1 = cuts_plus_border(column);
        x2 = cuts_plus_border(column+1);
        y1 = cuts_plus_border(row);
        y2 = cuts_plus_border(row+1);
        p = patch([x1 x1 x2 x2], [y1 y2 y2 y1], color);
        set(p, 'FaceAlpha', 0.75);
    end
end

