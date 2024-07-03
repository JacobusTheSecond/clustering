function [ nnidx_filtered, nndists_filtered ] = tw_filterNN( nnidx, nndists, keepPercent )
%TW_FILTERNN Summary of this function goes here
%   Detailed explanation goes here

    % ensure bounds
    keepPercent = min(keepPercent, 1);
    keepPercent = max(keepPercent, 0);

    dists = sort( nndists(~isinf(nndists)) );
    radius_idx = floor(numel(dists) * keepPercent);
    if radius_idx == 0
        nnidx_filtered = [];
        nndists_filtered = [];
    else
        filter_radius = dists(radius_idx);
        filtered_idx = nndists <= filter_radius;
        
        nnidx_filtered = nan(size(nnidx));
        nndists_filtered = inf(size(nndists));
        
        nnidx_filtered(filtered_idx) = nnidx(filtered_idx);
        nndists_filtered(filtered_idx) = nndists(filtered_idx);
    end
end

