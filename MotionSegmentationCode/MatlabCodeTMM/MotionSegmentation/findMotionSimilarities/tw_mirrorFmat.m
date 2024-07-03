function [ fmat_mirror ] = tw_mirrorFmat( fmat )
%TW_MIRRORFMAT Summary of this function goes here
%   Detailed explanation goes here

    fmat_mirror = zeros(size(fmat));
    
    feature_size = size(fmat, 1);
    offset_size = 5*3;
    offset_count = feature_size / offset_size;
    
    for i = 0:offset_count-1
        offset = i*offset_size;
        fmat_mirror(1+offset, :) = -fmat(4+offset, :);
        fmat_mirror(2+offset, :) = fmat(5+offset, :);
        fmat_mirror(3+offset, :) = fmat(6+offset, :);
        fmat_mirror(4+offset, :) = -fmat(1+offset, :);
        fmat_mirror(5+offset, :) = fmat(2+offset, :);
        fmat_mirror(6+offset, :) = fmat(3+offset, :);
        fmat_mirror(7+offset, :) = -fmat(7+offset, :);
        fmat_mirror(8+offset, :) = fmat(8+offset, :);
        fmat_mirror(9+offset, :) = fmat(9+offset, :);
        fmat_mirror(10+offset, :) = -fmat(13+offset, :);
        fmat_mirror(11+offset, :) = fmat(14+offset, :);
        fmat_mirror(12+offset, :) = fmat(15+offset, :);
        fmat_mirror(13+offset, :) = -fmat(10+offset, :);
        fmat_mirror(14+offset, :) = fmat(11+offset, :);
        fmat_mirror(15+offset, :) = fmat(12+offset, :);
    end
end

