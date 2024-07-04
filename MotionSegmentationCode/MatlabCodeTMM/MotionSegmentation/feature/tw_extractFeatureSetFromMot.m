function [ featureSet ] = tw_extractFeatureSetFromMot( skel, mot, feature_type, frame_offsets)
%TW_EXTRACTFEATURESETFROMMOT Summary of this function goes here
%   Detailed explanation goes here

    motMats    = convertMot2DBmat(skel, mot);

    featureSet = tw_extractFeatureSetFromDB(motMats, feature_type, frame_offsets);

    if any(isnan(featureSet(:)))
        error('Specified skel does not support specified feature set!');
    end
end

