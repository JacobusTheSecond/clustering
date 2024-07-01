function [ cuts ] = tw_findActivities( nnidx, nndists, options )
%TW_FINDACTIVITIES Summary of this function goes here
%   Detailed explanation goes here

    radius_activities = options.radius_activities;
    nn_min_blob_size = options.nn_min_blob_size;
    min_zero_deltas_cut = options.min_zero_deltas_cut;
    filter_diagonal_width = options.filter_diagonal_width;
    use_old_activity_search = options.use_old_activity_search;
    DEBUG_FIND_MOTION_ACTIVITIES = options.DEBUG_FIND_MOTION_ACTIVITIES;

    [nnidx_filtered nndists_filtered] = tw_filterRadius(radius_activities, nnidx, nndists);
    
    % get forward and backward cuts
    if (use_old_activity_search)
        [~,cutsf,cutsb] = tw_countNNperSquare_old(nnidx_filtered);
    else
        [cutsf, cutsb] = tw_countNNperSquare(nnidx_filtered, nn_min_blob_size, min_zero_deltas_cut, filter_diagonal_width);
    end
        
    % join forward and backward cuts
    cuts = sort([cutsf cutsb]);
    
    if (DEBUG_FIND_MOTION_ACTIVITIES)
        [ nnidx_filtered_diag, nndists_filtered_diag ] = tw_filterDiagonal( filter_diagonal_width, nnidx_filtered, nndists_filtered );
        frame_count = size(nnidx_filtered, 2);
        figure(21);
        tw_selfSimilarityPlot(size(nnidx_filtered_diag, 2), options.radius, nnidx_filtered_diag, nndists_filtered_diag);
        tw_plotCutsColor(cutsf, 'g');
        title('Forward Cuts');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(22);
        tw_selfSimilarityPlot(size(nnidx_filtered_diag, 2), options.radius, nnidx_filtered_diag, nndists_filtered_diag);
        tw_plotCutsColor(cutsf, 'g');
        tw_plotForwardScanlines(cutsf, frame_count, 'r');
        title('Forward Cuts Scanlines');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(23);
        tw_selfSimilarityPlot(size(nnidx_filtered_diag, 2), options.radius, nnidx_filtered_diag, nndists_filtered_diag);
        tw_plotCutsColor(cutsb, 'm');
        title('Backward Cuts');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(24);
        tw_selfSimilarityPlot(size(nnidx_filtered_diag, 2), options.radius, nnidx_filtered_diag, nndists_filtered_diag);
        tw_plotCutsColor(cutsb, 'm');
        tw_plotBackwardScanlines(cutsb, frame_count, 'r');
        title('Backward Cuts Scanlines');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(25);
        tw_selfSimilarityPlot(size(nnidx, 2), options.radius, nnidx, nndists);
        tw_plotCutsColor(cutsf, 'g');
        tw_plotCutsColor(cutsb, 'm');
        title('All Cuts');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(26);
        tw_selfSimilarityPlot(size(nnidx, 2), options.radius, nnidx, nndists);
        title('Sparse SSM');
        axis image;
        set(gcf,'position',[100 0 900 900]);
        xlabel('Frame Index');
        ylabel('Frame Index');
        figure(1);
    end
end
