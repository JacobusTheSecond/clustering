function [] = tw_plotForwardScanlines( cutsf, frame_count, color )
%TW_PLOTFORWARDSCANLINES Summary of this function goes here
%   Detailed explanation goes here

    cutsf = sort(cutsf);

    start_x = 1;
    cutsf_idx = 1;
    start_x_list = zeros(1, frame_count);
    end_x_list = zeros(1, frame_count);
    start_y_list = zeros(1, frame_count);
    end_y_list = zeros(1, frame_count);
    for i=1:frame_count
        start_x_list(i) = start_x;
        end_x_list(i) = i;
        start_y_list(i) = i;
        end_y_list(i) = i;
        
        if (cutsf(cutsf_idx) == i)
            start_x = i;
            cutsf_idx = min(numel(cutsf), cutsf_idx + 1);
        end
    end
    
    prune = 1:4:frame_count;
    
    start_x_list = start_x_list(prune);
    end_x_list = end_x_list(prune);
    start_y_list = start_y_list(prune);
    end_y_list = end_y_list(prune);
    
    line(vertcat(start_x_list, end_x_list), vertcat(start_y_list, end_y_list), 'Color', color, 'LineWidth', 0.5);
end