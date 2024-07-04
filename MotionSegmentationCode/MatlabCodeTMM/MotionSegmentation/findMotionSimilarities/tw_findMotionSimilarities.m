function [ nnidx, nndists, nnidx_mirror, nndists_mirror, time ] = tw_findMotionSimilarities( fmat, fmat_mirror, options )
%TW_FINDMOTIONSIMILARITIES Summary of this function goes here
%   Detailed explanation goes here
    time = 0;
    nnidx_mirror = [];
    nndists_mirror = [];

    frame_offsets = options.frame_offsets;
    use_normal_weights = options.use_normal_weights;
    normal_weights_sigma = options.normal_weights_sigma;
    k = options.k;
    radius = options.radius;
    use_mirror_motion = options.use_mirror_motion;
    
    
    %% apply normal distribution weights
    if use_normal_weights
        [ fmat ] = tw_applyNormalWeights(fmat, normal_weights_sigma, frame_offsets);
        
        if use_mirror_motion
            [ fmat_mirror ] = tw_applyNormalWeights(fmat_mirror, normal_weights_sigma, frame_offsets);
        end
    end
    
    %% get k nearest neighbors
    if use_mirror_motion
        [ nnidx, nndists, nnidx_mirror, nndists_mirror, time ] = tw_getMirrorNN( fmat, fmat_mirror, k, radius );
       else
        [ nnidx, nndists, time ] = tw_getNN( fmat, k, radius );
    end
    
end

