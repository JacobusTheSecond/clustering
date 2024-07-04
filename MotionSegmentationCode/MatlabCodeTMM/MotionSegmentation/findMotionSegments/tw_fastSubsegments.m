function [ all_paths, all_dists ] = tw_fastSubsegments(A_cut, sf_cut, min_path_length, nnidx_cut, nndists_cut)
%TW_FASTSUBSEGMENTS Summary of this function goes here
%   Detailed explanation goes here

    k_nn = size(nnidx_cut, 1);

    % only do anything if we have nn's
    A_strong_connected  = A_cut(2:end-1, 2:end-1);
    % Since A_strong_connected is a DAG, it is allowed to make A_strong_connected symmetric
    % by adding the transpose. No original element of A_strong_connected is changed.
    A_strong_connected = A_strong_connected + A_strong_connected';

    % Determine for each entry in the SSM its connected component
    comps = components_mex(A_strong_connected);
    comps_cut = reshape(comps, size(nnidx_cut));

    % Number of connected components == max(comps)
    % Number of elements inside a connected component
    comps_size = hist(comps, max(comps));
    
    % Indices of CC with more than one element
    % The indices equal the CC identifiers
    % nan components have size 1
    comps_valid = find(comps_size > 1);
    
    all_paths = cell( numel(comps_valid), 1 );
    all_dists = zeros( numel(comps_valid), 1 );
    result_idx = 1;

    for comp = comps_valid
        % get start frame and end frame
        % find all indexes for the entries belonging to the current CC
        [~, cols] = find( comps_cut == comp );
        sf_comp = cols(1);
        ef_comp = cols(end);

        % The CC is completely contained in a frame and therefore
        % irrelevant.
        if sf_comp == ef_comp
            continue
        end
                
        % compute start and end node
        s_dist = nndists_cut(:, sf_comp);
        %e_dist = nndists_cut(:, ef_comp);
        e_dist = ones(k_nn, 1) * 0.0001;
        
        s_dist(s_dist == 0) = 0.0001;
        %e_dist(e_dist == 0) = 0.0001;
        s_dist(comps_cut(:, sf_comp) ~= comp) = 0;
        e_dist(comps_cut(:, ef_comp) ~= comp) = 0;
        
        % clear old entries
        A_cut(1, :) = 0;
        A_cut(:, end) = 0;
        
        % set start and end node
        start_entry = ((sf_comp - 1) * k_nn + 1) + 1;
        end_entry = ((ef_comp - 1) * k_nn + 1) + 1;
        A_cut( 1, (start_entry:start_entry+k_nn-1) ) = s_dist;
        A_cut( (end_entry:end_entry+k_nn-1), end ) = e_dist;
        
        % Returns the distance (L) and the predecessor (pred) 
        % for each of the vertices along the shortest path 
        % from 1 to every other vertex in the graph.
        [L, pred] = dag_sp(A_cut, 1);

        length = 0;
        node = pred(end);

        % No path from the start node to the end node exists
        if node==0
            continue             
        end

        % compute path length
        % while node doesn't equal start node
        while node~=1
            node = pred(node);
            length = length + 1;
        end
        
        % The path isn't long enough, therefore the CC is too small to be
        % of interest.
        if length < min_path_length
            continue             
        end

        path = zeros(length, 2);

        node = pred(end);
        path_idx = 1;

        % while node doesn't equal start node
        while node~=1
            % node entry without start node
            node_entry = node - 1;
            
            %     node_entry = (node_cut_frame - 1) * k_nn + node_cut_frame_idx
            % <=> node_cut_frame = (node_entry - node_cut_frame_index) / k_nn + 1
            node_cut_frame = floor((node_entry - 1) / k_nn) + 1;
            node_frame = node_cut_frame + (sf_cut - 1);
            ref_node_frame = nnidx_cut( node_entry );

            path(path_idx, :) = [node_frame ref_node_frame];

            node = pred(node);
            path_idx = path_idx + 1;
        end

        all_paths(result_idx) = {path};
        all_dists(result_idx) = L(end);
        
        result_idx = result_idx + 1;
    end
    
    is_valid = all_dists > 0;
    all_paths = all_paths(is_valid);
    all_dists = all_dists(is_valid); 
end
