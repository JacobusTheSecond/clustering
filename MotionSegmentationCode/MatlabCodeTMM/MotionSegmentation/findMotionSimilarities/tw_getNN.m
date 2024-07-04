function [ nnidx, nndists,time ] = tw_getNN( fmat, k, radius )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % knn search
    frame_count = size(fmat, 2);
    handle_tree = ann_buildTree(fmat);
    k = min(k, frame_count);

    t0 = tic;
    t0 = toc;
    
    % nnidx: NN frame index for each frame
    % nndists: NN distance for each frame
    % The frames are compared to all other frames in the motion sequence
    [nnidx, nndists] = ann_queryTree(handle_tree, fmat, k,'search_sch', 'fr', 'radius', radius);
    
    nnidx = cast(nnidx, 'double');
    nnidx(nnidx==0) = nan;
    
    t1 = toc;
    
    ann_cleanup();
    
    time = t1-t0;
end

