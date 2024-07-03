function [allpaths, alldists] = tw_getSubsegments(nnidx, nndists, Afull, sf_global)
%TW_GETSUBSEGMENTS Summary of this function goes here
%   Detailed explanation goes here

    allpaths = {};
    alldists = [];

    k_nn = size(nnidx, 1);

    % only do anything if we have nn's
    if sum(sum(~isnan(nnidx))) > 2
        tmp  = Afull(2:end-1,2:end-1);
        % Since tmp is a DAG, it is allowed to make tmp symmetric
        % by adding the transpose. No original element of tmp is changed.
        tmp   = tmp+tmp';
        
        % Determine for each entry in the SSM its connected component
        comps = components_mex(tmp);
        
        
        compsmap = reshape(comps, size(nnidx));

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
            
            if sf~=ef
                % delete entries that do not belong to current component!
                Atmp = Afull;
                % filter out all entries that don't belong to 
                % the current CC
                mask = [curcomp;comps]~=curcomp;
                Atmp(mask,:) = 0;
                Atmp(:,mask) = 0;

                Ai = cutDAGMatrix(Atmp, k_nn, sf, ef, nndists(:, sf), ones(1, k_nn));

                % Returns the distance (L) and the predecessor (pred) 
                % for each of the vertices along the shortest path 
                % from 1 to every other vertex in the graph.  
                [L, pred]      = dag_sp(Ai, 1);
                
                % Returns paths sorted by distances:
                % Each path row consists of the frame index relative to the
                % cut out SSM part and the NN frame index of the path node.
                % The first row represents the end node and the last row 
                % the start node.
                % Returns also the distance from the start node to
                % the end node of each path.
                [paths, dists] = getTrellisPaths_acm(nnidx(:, sf:ef), L, pred);
                % Merge paths according to the index distance of the start
                % frames. A window size of 10 is used on both sides.
                % The path with smallest distance is used, all other paths
                % within the window are removed.
                [paths, dists] = deleteDoublePaths(paths, dists);

                for cp = 1:numel(paths)
                    % make frame index in paths absolute
                    paths{cp}(:, 1) = paths{cp}(:, 1) + sf - 1;
                end

                allpaths = [allpaths; paths];
                alldists = [alldists; dists];
            end      
        end

        % filter out empty paths
        allpaths = allpaths(~cellfun('isempty',allpaths));
        
        for c = 1:numel(allpaths)
            curpath = allpaths{c};
            curpath(:,1) = curpath(:,1) + sf_global - 1;
            allpaths{c} = curpath;
        end
    end
end