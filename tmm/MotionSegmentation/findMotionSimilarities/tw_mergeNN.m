function [ nnidx_merge, nndists_merge ] = tw_mergeNN( nnidx, nndists, nnidx_mirror, nndists_mirror, is_cut_to_k )
%TW_MERGENN Summary of this function goes here
%   Detailed explanation goes here
    
    [k_nn, frame_count] = size(nnidx);

    if (is_cut_to_k)
        k_merge = k_nn;
    else
        k_merge = k_nn * 2;
    end
    
    nnidx_merge = nan(k_merge, frame_count);
    nndists_merge = inf(k_merge, frame_count);

    for i=1:frame_count
        temp = [nnidx(:,i) nndists(:,i); nnidx_mirror(:,i) nndists_mirror(:,i)];
        % sort by index primary and by distance secondary
        temp_sorted = sortrows(temp, [1 2]);

        % remove dublicate entries
        dublicate_idx = [false; diff(temp_sorted(:,1))<(10^-6)];
        [temp_sorted, ~] = removerows(temp_sorted,'ind', dublicate_idx);

        % sort by distance primary and by index secondary
        temp_sorted = sortrows(temp_sorted, [2 1]);

        merge_range = 1:min( size(temp_sorted, 1), k_merge );
        nnidx_merge(merge_range, i) = temp_sorted(merge_range, 1);
        nndists_merge(merge_range, i) = temp_sorted(merge_range, 2);
    end
end

