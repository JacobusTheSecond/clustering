function [ nnidx_sym, nndists_sym ] = tw_makeSymmetric( nnidx, nndists )
%TW_MAKESYMMETRIC Summary of this function goes here
%   Detailed explanation goes here
    tic;
    fprintf('Make Symmetric. ')

    [k_nn, frame_count] = size(nnidx);
    nnidx_sym = nan(k_nn, frame_count);
    nndists_sym = inf(k_nn, frame_count);
   
    for i=1:frame_count
       is_entry = ~isnan( nnidx(:, i) );
       ref_idx = nnidx(is_entry, i);
       ref_dists = nndists(is_entry, i);
       
       ref_frames = nnidx(:, ref_idx);
       is_valid = logical(sum(ref_frames == i));
       valid_idx = ref_idx(is_valid);
       valid_dists = ref_dists(is_valid);
       valid_count = numel(valid_idx);
       
       nnidx_sym(1:valid_count, i) = valid_idx;
       nndists_sym(1:valid_count, i) = valid_dists;
    end
    
    toc;
    
end

