function [ nnidx_filtered, nndists_filtered ] = tw_filterRadius( radius, nnidx, nndists )
%TW_FILTERRADIUS Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_count] = size(nnidx);
    is_in_radius = nndists <= radius;
    nnidx_filtered = nan(k_nn, frame_count);
    nndists_filtered = inf(k_nn, frame_count);
    nnidx_filtered(is_in_radius) = nnidx(is_in_radius);
    nndists_filtered(is_in_radius) = nndists(is_in_radius);
end

