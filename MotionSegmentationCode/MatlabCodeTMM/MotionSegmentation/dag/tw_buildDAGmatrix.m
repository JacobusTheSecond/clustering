function A = tw_buildDAGmatrix(nnidx, nndists, allowedSteps, step_over_frame_penalty)
%TW_BUILDDAGMATRIX Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_count] = size(nnidx);
    allowed_steps_count = size(allowedSteps, 1);

    % Larger steps can connect components and are therefore included.
    % They should however not be used inside a connected component.
    % Therefore are penalty weight is added for each skipped frame.
    step_penalties = zeros(1, allowed_steps_count);
    for i = 1:allowed_steps_count
        vstep_penalty = max( 0, allowedSteps(i,1) - 1 ) * step_over_frame_penalty;
        hstep_penalty = max( 0, allowedSteps(i,2) - 1 ) * step_over_frame_penalty;

        step_penalties(i) = vstep_penalty + hstep_penalty;
    end

    % Set minimum distance to a non zero value.
    % We use a sparse matrix later, therefore zero means no connection.
    nndists( nndists == 0 ) = 0.01;

    % Prepare arrays for the sparse Adjacency Matrix
    entry_cols = frame_count * allowed_steps_count + 1; % start node
    from = zeros(k_nn, entry_cols);
    to = zeros(k_nn, entry_cols);
    weights = zeros(k_nn, entry_cols);
    
    % Construct graph beginning    
    % from the start node
    from(:, 1) = 1;
    % to the k NN of the first frame
    to(:, 1) = 2:k_nn+1;
    % use distance as weight
    s_dist = nndists(:, 1);
    s_dist(isinf(s_dist)) = 0;
    weights(:, 1) = s_dist;

    current_entry_col = 2;

    % for each frame
    for frame = 1:frame_count
        % for each step
        for step = 1:allowed_steps_count

            vstep = allowedSteps(step, 1);
            hstep = allowedSteps(step, 2);

            hframe = frame + hstep;
            % if hframe is still a frame
            if hframe > 0 && hframe <= frame_count
                % Number of nodes before frame nodes in the graph
                % Start Node + Preceding Nodes
                graph_nodes_before_frame = 1 + (frame - 1) * k_nn;
                % Number of nodes before hframe nodes in the graph
                % Start Node + Preceding Nodes
                graph_nodes_before_hframe = 1 + (hframe - 1) * k_nn;

                % frame NN + vstep
                frame_nnidx_vstep = nnidx(:, frame) + vstep;
                hframe_nnidx = nnidx(:, hframe);
                
                % nearest neighbors which are reachable with the current step
                [is_reachable, is_reachable_hframe_idx] = ismember(frame_nnidx_vstep, hframe_nnidx);
                % reachable nodes index in hframe
                reachable_hframe_idx = is_reachable_hframe_idx(is_reachable);
                
                penalty_weight = step_penalties(step);

                from(:, current_entry_col) = graph_nodes_before_frame + (1:k_nn);            
                to(is_reachable, current_entry_col) = graph_nodes_before_hframe + reachable_hframe_idx;
                weights(is_reachable, current_entry_col) = nndists(reachable_hframe_idx, hframe) + penalty_weight;

                current_entry_col = current_entry_col + 1;
            end
        end
    end

    % some entries in 'to' are zero and all corresponding 'weights' and
    % 'from' entries must be filtered out
    valid_entries = to > 0;
    from = from(valid_entries);
    to = to(valid_entries);
    weights = weights(valid_entries);

    % + start node + end node
    % Note that no connections to the end node has been set up.
    % The end node is inserted to allow efficient cutting and filling
    % in later steps.
    A_dim = frame_count * k_nn + 2;
    A = sparse(from, to, weights, A_dim, A_dim);
end