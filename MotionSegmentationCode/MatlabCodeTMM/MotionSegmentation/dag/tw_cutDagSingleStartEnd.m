function [ A_cut, nnidx_cut, nndists_cut ] = tw_cutDagSingleStartEnd( A, sf, ef, is_cut_to_range, nnidx, nndists )
%TW_CUTDAGSINGLESTARTEND Summary of this function goes here
%   Detailed explanation goes here

    k_nn = size(nnidx, 1);
    
    range = sf:ef;
    nnidx_cut = nnidx(:, range);
    nndists_cut = nndists(:, range);
    
    if (is_cut_to_range)
        smaller = nnidx_cut < sf;
        larger = nnidx_cut > ef;
        nnidx_cut(smaller) = nan;
        nnidx_cut(larger) = nan;
        nndists_cut(smaller) = inf;
        nndists_cut(larger) = inf;
    end

    s_dist = nndists_cut(:, 1);
    %e_dist = nndists_cut(:, end);
    e_dist = ones(1, k_nn) * 0.0001;
    
    s_dist(s_dist == 0) = 0.01;
    e_dist(e_dist == 0) = 0.01;
    s_dist(isinf(s_dist)) = 0;
    e_dist(isinf(e_dist)) = 0;
        
    s_idx = ((sf - 1) * k_nn + 1) + 1; % + start node
    e_idx = (ef * k_nn) + 1; % + start node
    
    % leave space for start and end node
    %entry_range = s_idx-1:e_idx+1;
    
    min_idx = s_idx;
    max_idx = e_idx;
    [row,col,v] = find(A);
    
    row_range_idx = (row >= min_idx & row <= max_idx);
    row = row(row_range_idx);
    col = col(row_range_idx);
    v = v(row_range_idx);
    
    col_range_idx = (col >= min_idx & col <= max_idx);
    row = row(col_range_idx);
    col = col(col_range_idx);
    v = v(col_range_idx);
    
    row = (row - min_idx + 1) + 1; % leave space for start node (index starts at 2)
    col = (col - min_idx + 1) + 1; % leave space for start node (index starts at 2)
    % leave space for start and end node
    count = ((max_idx - min_idx) + 1) + 1 + 1;
    A_cut = sparse(row, col, v, count, count, numel(v));
    
%     A_cut = A(entry_range, entry_range);
%     
%     % clear start and end node
%     A_cut(1, :) = 0;
%     A_cut(:, 1) = 0;
%     A_cut(end, :) = 0;
%     A_cut(:, end) = 0;
    
    % set start and end node
    A_cut(1, 2:k_nn+1) = s_dist;
    A_cut(end-k_nn:end-1, end) = e_dist;
    
    if (is_cut_to_range)
        mask = [false; isnan( nnidx_cut(:) ); false];

        A_cut(:, mask) = 0;
        A_cut(mask, :) = 0;
    end
end

