function A = tw_buildDagFast(nnidx, nndists, allowedSteps, step_over_frame_penalty, global_frame_count)
%TW_BUILDDAGMATRIX Summary of this function goes here
%   Detailed explanation goes here

    [k_nn, frame_count] = size(nnidx);
    allowed_steps_count = size(allowedSteps, 1);

    % Larger steps can connect components and are therefore included.
    % They should however not be used inside a connected component.
    % Therefore are penalty weight is added for each skipped frame.
    vstep_penalty = max( 0, allowedSteps(:,1) - 1 ) * step_over_frame_penalty;
    hstep_penalty = max( 0, allowedSteps(:,2) - 1 ) * step_over_frame_penalty;
    step_penalties = vstep_penalty + hstep_penalty;

    % Set minimum distance to a non zero value.
    % We use a sparse matrix later, therefore zero means no connection.
    nndists( nndists == 0 ) = 0.0001;

    to_rows_base = repmat(1:frame_count, k_nn, 1);
    to_rows_base = to_rows_base(:)';
    to_columns_base = nnidx(:)';
    to_columns_idx_base = repmat((1:k_nn)', 1, frame_count);
    to_columns_idx_base = to_columns_idx_base(:)';
    to_dists_base = nndists(:)';
    
    % remove invalid to_columns entries
    valid_to_columns_base = ~isnan(to_columns_base); 
    to_rows_base = to_rows_base(valid_to_columns_base);
    to_columns_base = to_columns_base(valid_to_columns_base);
    to_columns_idx_base = to_columns_idx_base(valid_to_columns_base);
    to_dists_base = to_dists_base(valid_to_columns_base);
    
    % create check distance and check column index matrix
    to_check_index = tw_toLinearIndexRowWise(to_rows_base, to_columns_base, global_frame_count);
    to_check_index_ones = ones(1, numel(to_check_index));
    A_check_dist = sparse(to_check_index_ones, to_check_index, to_dists_base, 1, frame_count * global_frame_count, numel(to_check_index));
    A_check_columns_idx = sparse(to_check_index_ones, to_check_index, to_columns_idx_base, 1, frame_count * global_frame_count, numel(to_check_index));
    clear to_check_index;
    clear to_check_index_ones;
    
    % create from_rows and to_rows
    to_rows = repmat(to_rows_base, allowed_steps_count, 1);
    v_steps = repmat(allowedSteps(:,1), 1, numel(to_rows_base));
    from_rows = to_rows - v_steps;
    from_rows = from_rows(:)';
    to_rows = to_rows(:)'; 
    clear to_rows_base;
    clear v_steps;
    
    % create from_columns and to_columns
    to_columns = repmat(to_columns_base, allowed_steps_count, 1);
    to_columns_idx = repmat(to_columns_idx_base, allowed_steps_count, 1);
    h_steps = repmat(allowedSteps(:,2), 1, numel(to_columns_base));
    from_columns = to_columns - h_steps;
    from_columns = from_columns(:)';
    to_columns = to_columns(:)';
    to_columns_idx = to_columns_idx(:)';
    clear to_columns_base;
    clear to_columns_idx_base;
    clear h_steps;
    
    % create to_dists
    to_dists = repmat(to_dists_base, allowed_steps_count, 1);
    to_step_penalties = repmat(step_penalties, 1, numel(to_dists_base));
    to_dists = to_dists + to_step_penalties;
    to_dists = to_dists(:)';  
    clear to_dists_base;
    clear to_step_penalties;
    
    % remove invalid from_rows entries
    valid_from_rows = from_rows > 0;
    from_rows = from_rows(valid_from_rows);
    to_rows = to_rows(valid_from_rows);
    from_columns = from_columns(valid_from_rows);
    to_columns = to_columns(valid_from_rows);
    to_columns_idx = to_columns_idx(valid_from_rows);
    to_dists = to_dists(valid_from_rows);
    clear valid_from_rows;
    
    % remove invalid from_columns entries
    valid_from_columns = from_columns > 0;
    from_rows = from_rows(valid_from_columns);
    to_rows = to_rows(valid_from_columns);
    from_columns = from_columns(valid_from_columns);
    to_columns = to_columns(valid_from_columns);
    to_columns_idx = to_columns_idx(valid_from_columns);
    to_dists = to_dists(valid_from_columns);
    clear valid_from_columns;
    
    % remove invalid from_dists entries
    from_check_dist_index = tw_toLinearIndexRowWise(from_rows, from_columns, global_frame_count);
    from_dists = full( A_check_dist(from_check_dist_index) );
    valid_from_dists = from_dists ~= 0;
    from_rows = from_rows(valid_from_dists);
    to_rows = to_rows(valid_from_dists);
    from_columns = from_columns(valid_from_dists);
    to_columns = to_columns(valid_from_dists);
    to_columns_idx = to_columns_idx(valid_from_dists);
    to_dists = to_dists(valid_from_dists);
    clear from_check_dist_index;
    clear A_check_edge;
    clear valid_from_dists;
    
    % get from_columns_idx
    from_check_columns_idx_index = tw_toLinearIndexRowWise(from_rows, from_columns, global_frame_count);
    from_columns_idx = full( A_check_columns_idx(from_check_columns_idx_index) );
    clear from_check_columns_idx_index;
    clear A_check_edge;
    
    start_node = 1;
    end_node = 1;
    from = tw_toLinearIndexRowWise(from_rows, from_columns_idx, k_nn) + start_node;
    to = tw_toLinearIndexRowWise(to_rows, to_columns_idx, k_nn) + start_node;
    
    % + start node + end node
    % Note that no connections to the start node and to the end node have 
    % been set up.
    % The start node and end node are inserted to allow efficient cutting 
    % and filling in later steps.
    A_dim = frame_count * k_nn + start_node + end_node;
    A = sparse(from, to, to_dists, A_dim, A_dim);
end