function [ paths_valid, dists_valid, is_valid_symmetry ] = tw_checkSoftSymmetry( paths, dists, symmetry_coverage )
%TW_CHECKSOFTSYMMETRY Summary of this function goes here
%   Detailed explanation goes here

    is_valid_symmetry = false(numel(paths), 1);
    
    for i = 1:numel(paths)
        path = paths{i};
             
        for j = 1:numel(paths)
            % Do not have to filter out i == j,
            % since the path is compared to its mirrored version.
            
            mirrored_path = paths{j};
            
            [ frame_coverage, ref_frame_coverage ] = tw_getMirroredPathCoverage(path, mirrored_path);
            
            if frame_coverage >= symmetry_coverage && ref_frame_coverage >= symmetry_coverage;
                is_valid_symmetry(i) = true;
                break
            end
        end
    end

    paths_valid = paths(is_valid_symmetry);
    dists_valid = dists(is_valid_symmetry);
end

function [ frame_coverage, ref_frame_coverage ] = tw_getMirroredPathCoverage(path, mirrored_path)
    length = size(path, 1);
    
    frame_coverage = sum( ismember(path(:,1), mirrored_path(:,2)) ) / length;
    ref_frame_coverage = sum( ismember(path(:,2), mirrored_path(:,1)) ) / length;
end
