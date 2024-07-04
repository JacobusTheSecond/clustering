function [ skel, mot, name ] = tw_get_mot_by_index( index, base_path, skel_file, file_prefix, file_suffix )
%TW_GET_MOT_BY_INDEX Summary of this function goes here
%   Detailed explanation goes here
    skel_path = strcat(base_path, skel_file);
    file_path = sprintf('%s%s%02d%s', base_path, file_prefix, index, file_suffix);
    name = sprintf('%s%02d', file_prefix, index);
    
    
    [skel,mot] = readMocap(skel_path, file_path);
end

