function [ subcuts ] = tw_mergeMirrorCuts( subcuts_main, subcuts_mirror, min_cut_distance )
%TW_MERGEMIRRORCUTS Summary of this function goes here
%   Detailed explanation goes here
    
    if numel(subcuts_main) == 0
        subcuts = subcuts_mirror;
        return
    end

    subcuts = subcuts_main;
    subcuts_main_idx = 1;
    subcuts_main_count = numel(subcuts_main);
    
    lesser_neighbor = subcuts_main(subcuts_main_idx);
    greater_neighbor = subcuts_main(min(subcuts_main_idx+1, subcuts_main_count));
    
    for i = 1:numel(subcuts_mirror)
        subcut_mirror = subcuts_mirror(i);
        while (subcuts_main_idx < subcuts_main_count) ...
                && (greater_neighbor < subcut_mirror)
            subcuts_main_idx = subcuts_main_idx + 1;
            lesser_neighbor = subcuts_main(subcuts_main_idx);
            greater_neighbor = subcuts_main(min(subcuts_main_idx+1, subcuts_main_count));
        end
        
        if (abs(subcut_mirror - lesser_neighbor) > min_cut_distance) ...
                && (abs(subcut_mirror - greater_neighbor) > min_cut_distance)
            subcuts = [subcuts subcut_mirror];
        end
    end
end

