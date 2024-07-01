function [] = tw_selfSimilarityPlot(frame_count, radius, nnidx, nndists)
    %% construct self similarity matrix
    immat = nan(frame_count);
    for i=1:frame_count
        idx = ~isnan(nnidx(:, i));
        immat(i, nnidx(idx, i)) = nndists(idx, i);
    end
    
    %% create color map
    idx_range = 128;
    nan_color = [1 1 1];
    
    gray_map = vertcat(flipud(winter(idx_range)), nan_color);
    
    %% create image
    img_idx = floor(immat / radius * idx_range);
    img_idx(isnan(img_idx)) = idx_range + 1;
    img = ind2rgb(img_idx, gray_map);
    
    %% plot image
    image(img);
end