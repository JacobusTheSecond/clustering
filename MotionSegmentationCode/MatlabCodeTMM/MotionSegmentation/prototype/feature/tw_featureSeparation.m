function [ fmat_foot_l, fmat_foot_r, fmat_head, fmat_hand_l, fmat_hand_r ] = tw_featureSeparation( fmat, sf, ef )
%TW_FEATURESEPARATION Summary of this function goes here
%   Detailed explanation goes here

% id 4 name ltibia
% id 9 name rtibia
% id 17 name head
% id 21 name lwrist
% id 28 name rwrist
%[4,9,17,21,28]

    feature_count = size(fmat, 1) / 5;
    range = sf:ef;
    range_count = numel(range);
    fmat_foot_l = zeros(feature_count, range_count);
    fmat_foot_r = zeros(feature_count, range_count);
    fmat_head = zeros(feature_count, range_count);
    fmat_hand_l = zeros(feature_count, range_count);
    fmat_hand_r = zeros(feature_count, range_count);
    
    features = 1:3;    
    offset = [0 15 30];
    fmat_foot_l(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);
    
    offset = offset + numel(features);
    fmat_foot_r(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);

    offset = offset + numel(features);
    fmat_head(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);

    offset = offset + numel(features);
    fmat_hand_l(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);
    
    offset = offset + numel(features);
    fmat_hand_r(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);
end

