function [ fmat_foot, fmat_head, fmat_hand ] = tw_featureSeparation3Parts( fmat, sf, ef )
%TW_FEATURESEPARATION3PARTS Summary of this function goes here
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
    fmat_foot = zeros(feature_count * 2, range_count);
    fmat_head = zeros(feature_count, range_count);
    fmat_hand = zeros(feature_count * 2, range_count);
    
    offset = [0 15 30];
    features = 1:6;
    fmat_foot(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);
    
    offset = offset + numel(features);
    features = 1:3;
    fmat_head(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);

    offset = offset + numel(features);
    features = 1:6;
    fmat_hand(:, range) = fmat([offset(1) + features, offset(2) + features, offset(3) + features], range);

end
