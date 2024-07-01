function [ dtw_index ] = tw_toLinearIndexRowWise( row, column, column_size )
%TW_TOLINEARINDEXROWWISE Summary of this function goes here
%   Detailed explanation goes here

    dtw_index = (row - 1) * column_size + column;
end

