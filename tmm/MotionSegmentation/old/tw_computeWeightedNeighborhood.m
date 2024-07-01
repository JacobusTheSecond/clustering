function [w_idx w_dists immat] = tw_computeWeightedNeighborhood(fmat, idx, dists, minwidth, maxwidth)
%TW_COMPUTEWEIGHTEDNEIGHBORHOOD Summary of this function goes here
%   Detailed explanation goes here

frame_count = size(fmat,2);

w_idx = idx;
w_dists = dists;

posdiff = diff5point2(fmat,1);

distspf = normOfColumns(posdiff);
distspf = filterTimeline(distspf,11);

v_min = 0;
v_max = 40;

m = (maxwidth-minwidth)/(v_max-v_min);
m=2*m;

for i = 1:frame_count
    rcur = m*distspf(i) + minwidth;

    w_idx(w_dists(:,i) > rcur, i) = nan;
    w_dists(w_dists(:,i) > rcur, i) = inf;
end

idxperrow = sum(~isnan(w_idx), 2);
w_idx = w_idx(idxperrow > 0, :);
w_dists = w_dists(idxperrow > 0, :);


immat = nan(frame_count);

for i=1:frame_count
    w_idx_values_col_i = ~isnan( w_idx(:,i) );
    immat(i, w_idx(w_idx_values_col_i, i)) = dists(w_idx_values_col_i, i) ;
end

end