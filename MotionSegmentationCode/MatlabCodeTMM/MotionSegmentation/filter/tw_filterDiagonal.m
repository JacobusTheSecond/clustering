function [ nnidx_filtered, nndists_filtered ] = tw_filterDiagonal( diag_width, nnidx, nndists )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    frame_count = size(nnidx, 2);
    nnidx_filtered = nnidx;
    
    if nargin == 3
        nndists_filtered = nndists;
        for i = 1:frame_count
            idx = abs(nnidx_filtered(:,i)-i) < diag_width;
            nnidx_filtered(idx, i) = nan;
            nndists_filtered(idx, i) = inf;
        end
    else
        for i = 1:frame_count
            idx = abs(nnidx_filtered(:,i)-i) < diag_width;
            nnidx_filtered(idx, i) = nan;
        end
    end
end