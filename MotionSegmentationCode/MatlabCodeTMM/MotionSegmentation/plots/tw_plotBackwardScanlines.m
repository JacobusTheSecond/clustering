function [ output_args ] = tw_plotBackwardScanlines( cutsb, frame_count, color )
%TW_PLOTBACKWARDSCANLINES Summary of this function goes here
%   Detailed explanation goes here

    cutsb = sort(cutsb);

    end_x = frame_count;
    cutsb_idx = numel(cutsb);
    start_x_list = zeros(1, frame_count);
    end_x_list = zeros(1, frame_count);
    start_y_list = zeros(1, frame_count);
    end_y_list = zeros(1, frame_count);
    for i=frame_count:-1:1
        start_x_list(i) = i;
        end_x_list(i) = end_x;
        start_y_list(i) = i;
        end_y_list(i) = i;
        if ~isempty(cutsb)
            if (cutsb(cutsb_idx) == i)
                end_x = i;
                cutsb_idx = max(1, cutsb_idx - 1);
            end
        end
    end
    
    prune = 1:4:frame_count;
    
    start_x_list = start_x_list(prune);
    end_x_list = end_x_list(prune);
    start_y_list = start_y_list(prune);
    end_y_list = end_y_list(prune);
    
    line(vertcat(start_x_list, end_x_list), vertcat(start_y_list, end_y_list), 'Color', color, 'LineWidth', 0.5);
end

