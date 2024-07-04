function [] = tw_mirrorPlot(frame_count, radius, nnidx, nndists, nnidx_mirror, nndists_mirror)
%TW_MIRRORPLOT Summary of this function goes here
%   Detailed explanation goes here

    %% construct self similarity matrices
    immat = nan(frame_count);
    for i=1:frame_count
        idx = ~isnan(nnidx(:, i));
        immat(i, nnidx(idx, i)) = nndists(idx, i);
    end
    
    immat_mirror = nan(frame_count);
    for i=1:frame_count
        idx = ~isnan(nnidx_mirror(:, i));
        immat_mirror(i, nnidx_mirror(idx, i)) = nndists_mirror(idx, i);
    end
    
    immat_merge = nan(frame_count);
    for i=1:frame_count
        idx = ~isnan(nnidx(:, i));
        immat_merge(i, nnidx(idx, i)) = nndists(idx, i);
        idx = ~isnan(nnidx_mirror(:, i));
        immat_merge(i, nnidx_mirror(idx, i)) = nndists_mirror(idx, i);
    end
    
    %% create color maps
    idx_range = 128;
    nan_color = [ 1 1 1 ];
    
    gray_map = vertcat(flipud(winter(idx_range)), nan_color);
    mirror_map = vertcat(flipud(autumn(idx_range)), nan_color);
        
    %% create merge image
    merge_idx = floor(immat_merge / radius * idx_range);
    merge_idx(isnan(merge_idx)) = idx_range + 1;
    merge_img = ind2rgb(merge_idx, gray_map);
    
    %% extract mirrored part in merge image
    diffmat_mirror = immat_merge - immat_mirror;
    immat_mirror_part = nan(size(immat_merge));
    immat_mirror_part((diffmat_mirror == 0)) = immat_merge((diffmat_mirror == 0));
        
    mirror_part_idx = floor(immat_mirror_part / radius * idx_range);
    mirror_nans = isnan(mirror_part_idx);
    mirror_part_idx(mirror_nans) = idx_range + 1;
    mirror_part_img = ind2rgb(mirror_part_idx, mirror_map);
    
    %% compose final image
    mirror_values_3 = false(size(mirror_nans, 1), size(mirror_nans, 2), 3);
    mirror_values = ~mirror_nans;
    mirror_values_3(:,:,1) = mirror_values;
    mirror_values_3(:,:,2) = mirror_values;
    mirror_values_3(:,:,3) = mirror_values;
    final_img = merge_img;
    final_img(mirror_values_3) = mirror_part_img(mirror_values_3);
    
    %% plot final image
    image(final_img);
end