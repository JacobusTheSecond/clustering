function [ A_cut, nnidx_cut, nndists_cut ] = tw_createCutDagSingleStartEnd( sf, ef, is_cut_to_range, nnidx, nndists, allowedSteps, step_over_frame_penalty)
%TW_CUTDAGSINGLESTARTEND Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_count] = size(nnidx);
    
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
        
    A_cut = tw_buildDagFast(nnidx_cut, nndists_cut, allowedSteps, step_over_frame_penalty, frame_count);

    s_dist = nndists_cut(:, 1);
    %e_dist = nndists_cut(:, end);
    e_dist = ones(1, k_nn) * 0.0001;

    s_dist(s_dist == 0) = 0.01;
    e_dist(e_dist == 0) = 0.01;
    s_dist(isinf(s_dist)) = 0;
    e_dist(isinf(e_dist)) = 0;
    
    % set start and end node
    A_cut(1, 2:k_nn+1) = s_dist;
    A_cut(end-k_nn:end-1, end) = e_dist;
end

