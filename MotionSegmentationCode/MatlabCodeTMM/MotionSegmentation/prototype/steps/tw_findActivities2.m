function [ cuts ] = tw_findActivities2( nnidx, nndists, radius_s1, sampling_rate )
%TW_FINDACTIVITIES2 Summary of this function goes here
%   Detailed explanation goes here

    tic   
    %[nnidx_filtered ~] = tw_filterRadius(radius_s1, nnidx, nndists);
    nnidx_filtered = nnidx;
    
    % get forward and backward cuts
    [cutsf, cutsb] = tw_countNNperSquare2(nnidx_filtered, sampling_rate);

    % join forward and backward cuts
    cuts = sort([cutsf cutsb]);
    
    fprintf('Find Activities completed. ')
    toc
end