function [ fmat, fmat_mirror ] = tw_createFeatureSet( skel, mot, options )
%TW_CREATEFEATURESET Summary of this function goes here
%   Detailed explanation goes here

    feature_set = options.feature_set;
    frame_offsets = options.frame_offsets;
    use_mirror_motion = options.use_mirror_motion;
    use_feature_projection = options.use_feature_projection;
    feature_projection_k = options.feature_projection_k;
    
    fmat_mirror = [];

    %% extract feature set
    if use_feature_projection
        options.datatype = 'mocapProjected';
        fprintf('Starting Kernel based projection.\n')

%        tic;
        motn   = prepareMotForQuery(skel, mot);
        fmat   = tw_extractFeatureSetFromMot(skel, motn, feature_set, frame_offsets);
        projection_options.VIZPROJECTION = true;
        projection_options.k = feature_projection_k;
        [pfmat,meta] = computeProjectionKernelBased(fmat, projection_options);

%         timings.projection = toc;
%         timings.projectionMeta = meta;

%         mot.data = pfmat;
%         fmat = mot.data;
        fmat = pfmat;
        
        if use_mirror_motion
%             [skel_mirror, mot_mirror] = mirrorMot(skel, mot);
%             motn_mirror   = prepareMotForQuery(skel_mirror, mot_mirror);
%             fmat_mirror   = tw_extractFeatureSetFromMot(skel_mirror, motn_mirror, feature_set, frame_offsets);
%             [fmat_mirror, meta] = computeProjectionKernelBased(fmat_mirror);
            
            [ fmat_mirror ] = tw_mirrorFmat( fmat );
        end

        fprintf('Finished projecting points!\n')
    else
        fmat = tw_extractFeatureSetFromMot(skel, mot, feature_set, frame_offsets);
        
        if use_mirror_motion
            [skel_mirror, mot_mirror] = mirrorMot(skel, mot);
            fmat_mirror = tw_extractFeatureSetFromMot(skel_mirror, mot_mirror, feature_set, frame_offsets);
        end
    end
end

