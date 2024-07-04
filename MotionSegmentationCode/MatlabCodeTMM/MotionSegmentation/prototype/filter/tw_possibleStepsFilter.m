function [ nnidx_filtered, nndists_filtered ] = tw_possibleStepsFilter(cuts, nnidx, nndists)
%TW_POSSIBLESTEPSFILTER Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_count] = size(nnidx);
    
    is_forward = false(k_nn, frame_count);
    is_backward = false(k_nn, frame_count);

    frames2skip_begin = 0;
    allowedSteps = [
                     1 1; ...
%                      1 0; ...
%                      0 1; ...
                     1 2; ...
                     2 1; ...
                     2 2; ...
%                      5 5; ...
%                      10 10; ...
                   ];
               
    Afull = tw_buildDAGmatrix(nnidx, nndists, frames2skip_begin, allowedSteps);
    
    % for all motion primitives
    for i=1:numel(cuts)+1
        if i==1
           sf = 1;
        else
           sf = cuts(i-1) + 1;
        end
        if i==numel(cuts)+1
           ef = frame_count;
        else
           ef = cuts(i);
        end

        [A_strip, strip_idx, strip_dists] = cutDAGMatrixRect(Afull, sf, ef, 1, frame_count, k_nn, nnidx, nndists);
        %[A_strip, strip_idx, strip_dists] = cutDAGMatrixSquare(Afull, sf, ef, k_nn, nnidx, nndists);
        [is_forward_strip, is_backward_strip] = tw_isReachable(frame_count, sf, A_strip, strip_idx, strip_dists);
    
        is_forward(:, sf:ef) = is_forward_strip;
        is_backward(:, sf:ef) = is_backward_strip;
    end
    
    nnidx_filtered = nan(k_nn, frame_count);
    nndists_filtered = inf(k_nn, frame_count);
    
    nnidx_filtered(is_forward) = nnidx(is_forward);
    %nnidx_filtered(is_backward) = nnidx(is_backward);
    
    nndists_filtered(is_forward) = nndists(is_forward);
    %nndists_filtered(is_backward) = nndists(is_backward);
end

function [is_forward, is_backward] = tw_isReachable(frame_count, offset, A_strip, strip_idx, strip_dists)
%TW_GETSUBSEGMENTS Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_range] = size(strip_idx);
    
    is_forward = false(1, k_nn * frame_range);
    is_backward = false(1, k_nn * frame_range);

    % only do anything if we have nn's
    if sum(sum(~isnan(strip_idx))) > 2
        tmp  = A_strip(2:end,2:end);
        % Since tmp is a DAG, it is allowed to make tmp symmetric
        % by adding the transpose. No original element of tmp is changed.
        tmp   = tmp+tmp';
        
        % Determine for each entry in the SSM its connected component
        comps = components_mex(tmp);
        
        compsmap = reshape(comps, size(strip_idx));

        % Number of connected components == max(comps)
        % Number of elements inside a connected component
        compsize = hist(comps, max(comps));
        % Indices of CC with more than one element
        % The indices equal the CC identifiers 
        relcomponents = find(compsize>1);

        % for each connected component with more than one element
        for curcomp = relcomponents

            % get start frame and end frame
            % find all indexes for the entries belonging to the current CC
            [rows,cols] = find(compsmap==curcomp);
            sf = cols(1);
            ef = cols(end);
            
            %TODO remove
%             t_f = find(compsmap==curcomp);
%             t_idx = nan(k_nn, frame_range);
%             t_dists = inf(k_nn, frame_range);
%             t_idx(t_f) = strip_idx(t_f);
%             t_dists(t_f) = strip_dists(t_f);
%             tf_idx = nan(k_nn, frame_count);
%             tf_dists = inf(k_nn, frame_count);
%             tf_idx(:, offset:(offset+frame_range-1)) = t_idx;
%             tf_dists(:, offset:(offset+frame_range-1)) = t_dists;
%             figure(3);
%             tw_selfSimilarityPlot(frame_count, 64, tf_idx, tf_dists);
%             figure(1);
%             
%             dbstop in title
%             title('abc');
            
            if sf~=ef
                % delete entries that do not belong to current component!
                Atmp = A_strip;
                % filter out all entries that don't belong to 
                % the current CC
                mask = [curcomp;comps] ~= curcomp;
                Atmp(mask,:) = 0;
                Atmp(:,mask) = 0;
                start_node_dists = strip_dists(:, sf);
                start_node_dists( compsmap(:, sf) ~= curcomp ) = inf;

                % forward reachable
                Ai = cutDAGMatrix(Atmp, k_nn, sf, ef, start_node_dists, ones(1, k_nn));
                Ai(isinf(Ai)) = 0;
                [d dt ft pred] = dfs(Ai, 1);
                pred = pred(2:end);
                forward_reach = find(pred);
                
                % backward reachable
                % use Atmp transpose
%                 Ai = cutDAGMatrix(Atmp', k_nn, sf, ef, strip_dists(:, sf), ones(1, k_nn));
%                 [d dt ft pred] = dfs(Ai, 1);
%                 pred = pred(2:end);
%                 % TODO flip pred
%                 backward_reach = find(pred);
                
                
                forward_reach = forward_reach + (sf-1) * k_nn;
%                 backward_reach = backward_reach + (sf-1) * k_nn;
                
                %TODO remove
                %t_f = find(compsmap==curcomp);
%                 t_f = forward_reach;
%                 t_idx = nan(k_nn, frame_range);
%                 t_dists = inf(k_nn, frame_range);
%                 t_idx(t_f) = strip_idx(t_f);
%                 t_dists(t_f) = strip_dists(t_f);
%                 tf_idx = nan(k_nn, frame_count);
%                 tf_dists = inf(k_nn, frame_count);
%                 tf_idx(:, offset:(offset+frame_range-1)) = t_idx;
%                 tf_dists(:, offset:(offset+frame_range-1)) = t_dists;
%                 figure(3);
%                 tw_selfSimilarityPlot(frame_count, 64, tf_idx, tf_dists);
%                 figure(1);
% 
%                 dbstop in title
%                 title('abc');
                
                is_forward(forward_reach) = true;
%                 is_backward(backward_reach) = true;
            end      
        end
    end
    
    is_forward = reshape(is_forward, size(strip_idx));
    is_backward = reshape(is_backward, size(strip_idx));
end

