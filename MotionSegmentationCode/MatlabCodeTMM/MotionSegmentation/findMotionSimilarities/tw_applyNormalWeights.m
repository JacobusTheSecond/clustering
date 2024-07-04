function [ mfmat ] = tw_applyNormalWeights( fmat, sigma, offsets )
%TW_APPLYNORMALWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

    mfmat = zeros(size(fmat));
    
    weights = normpdf(offsets, 0, sigma);
    
    % the zero offset should have weight one
    %weights = weights * 1/weights(offsets == 0);
    
    scale = numel(offsets) / sum(weights);
    weights = scale * weights;
    weights = sqrt(weights);
    fsize = size(fmat, 1) / numel(offsets);        
    for i = 1:numel(offsets)
        start_index = 1 + fsize * (i-1);
        end_index = fsize * i;

        mfmat(start_index:end_index, :) = fmat(start_index:end_index, :) * weights(i);
    end
end

