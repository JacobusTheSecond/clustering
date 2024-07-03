function [ nnc, nncb, cuts, cutsb ] = tw_countNNperSquare3( nnidx, nndists, radius, sampling_rate )
%TW_COUNTNNPERSQUARE3 Summary of this function goes here
%   Detailed explanation goes here

[k_nn, frame_count] = size(nnidx);

potential_cut_forward = zeros(1, frame_count);

data_count = 5;
nnc = nan(data_count, frame_count);
nncb = nan(data_count, frame_count);

cuts = nan(1, frame_count);
cutsb = nan(1, frame_count);

[sort_idx, sort_permutation] = sort(nnidx);
sort_dists = nndists(sort_permutation);

% start frame
sf = 1;

%% find cuts forward 
for frame = 1:frame_count
    frame_idx = sort_idx(:, frame);
    frame_dists = sort_dists(:, frame);
    
    if frame==350 || frame==600 || frame==875
        figure(9);
        plot(frame_idx, frame_dists, '*');
        figure(1);

        dbstop in title
        title('abc')
    end
    
    [frame_pos] = find(frame == frame_idx);
    
    
    
    % remove neighbors that are smaller than sf and greater than i
%     smaller = frame_idx < sf;
%     greater = frame_idx > frame;
%     frame_idx(smaller)= nan;
%     frame_idx(greater) = nan;
%     frame_dists(smaller) = inf;
%     frame_dists(greater) = inf;
    
    sort(frame_idx);

    dists = frame_dists( ~isinf(frame_dists) );
    
    for d = 1:data_count
       nnc(d, frame) = sum(dists < (radius * d / data_count));
    end
end

%% find cuts backwards

% end frame
ef = frame_count;

for frame = frame_count:-1:1
    frame_idx = nnidx(:, frame);
    frame_dists = nndists(:, frame);

    % remove neighbors that are smaller than sf and greater than i
    smaller = frame_idx < frame;
    greater = frame_idx > ef;
    frame_idx(smaller)= nan;
    frame_idx(greater) = nan;
    frame_dists(smaller) = inf;
    frame_dists(greater) = inf;
    
    dists = frame_dists( ~isinf(frame_dists) );
    
    for d = 1:data_count
       nncb(d, frame) = sum(dists < (radius * d / data_count));
    end
    
    %nncb(frame) = nn_count_frame;
end

cuts = cuts(~isnan(cuts));
cutsb = cutsb(~isnan(cutsb));

end