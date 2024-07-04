function [ nnc, nncb, cuts, cutsb ] = tw_countNNperSquare2( nnidx, sampling_rate )
%TW_COUNTNNPERSQUARE2 Summary of this function goes here
%   Detailed explanation goes here

[k_nn, frame_count] = size(nnidx);

nnc = nan(1, frame_count);
nncb = nan(1, frame_count);

cuts = nan(1, frame_count);
cutsb = nan(1, frame_count);

% start frame
sf = 1;

nn_min_count = 50;
nn_min_max_diff = 50;
nn_max_min_diff = 50;
same_cut_candidate_max = 3;

activity_min = k_nn;
activity_max = 0;
activity_max_frame = nan;
has_peak = false;
same_cut_candidate_count = 0;
cut_candidate_min = k_nn;
cut_candidate_frame = nan;

%% find cuts forward 
for frame = 1:frame_count
    nn_frame = nnidx(:, frame);

    % remove neighbors that are smaller than sf and greater than i
    nn_frame(nn_frame < sf)= nan;
    nn_frame(nn_frame > frame) = nan;

    nn_count_frame = sum(~isnan(nn_frame));
    
    nnc(frame) = nn_count_frame;
    
    if has_peak
        if nn_count_frame == cut_candidate_min
            same_cut_candidate_count = same_cut_candidate_count + 1;
            if same_cut_candidate_count >= same_cut_candidate_max
                cuts(frame) = activity_max_frame;%cut_candidate_frame;
                sf = cut_candidate_frame;
                %activity_min = nn_count_frame;
                activity_max = nn_count_frame;
                activity_max_frame = nan;
                has_peak = false;
                same_cut_candidate_count = 0;
                cut_candidate_min = k_nn;
                cut_candidate_frame = nan;
            end
            
        elseif nn_count_frame < cut_candidate_min
            cut_candidate_min = nn_count_frame;
            cut_candidate_frame = frame;
            same_cut_candidate_count = 1;
            
        elseif (nn_count_frame - cut_candidate_min) > nn_min_max_diff
            cuts(frame) = activity_max_frame;%cut_candidate_frame;
            sf = cut_candidate_frame;
            %activity_min = cut_candidate_min;
            activity_max = nn_count_frame;
            activity_max_frame = nan;
            has_peak = false;
            same_cut_candidate_count = 0;
            cut_candidate_min = k_nn;
            cut_candidate_frame = nan;
        end
        
    elseif activity_max > nn_min_count && (activity_max - nn_count_frame) >  nn_max_min_diff
        has_peak = true;
        cut_candidate_min = nn_count_frame;
        cut_candidate_frame = frame;
        same_cut_candidate_count = 1;
        
    else
%         if (nn_count_frame < activity_min)
%             activity_min = nn_count_frame;
%         end
        if (nn_count_frame > activity_max)
            activity_max = nn_count_frame;
            activity_max_frame = frame;
        end
    end
end

%% find cuts backwards

% end frame
ef = frame_count;
nn_blob_size = 0;
zero_delta_count=0;

for frame = frame_count:-1:1
    nn_frame = nnidx(:, frame);

    % remove neighbors that are smaller than i and greater than ef
    nn_frame(nn_frame < frame) = nan;
    nn_frame(nn_frame > ef) = nan;
    
    nn_count_frame = sum(~isnan(nn_frame));
    nncb(frame) = nn_count_frame;
    % sum up all nearest neighbors from ef to current frame
    nn_blob_size = nn_blob_size + nn_count_frame;
 
    % if frame is a cut then ignore the current blob
    if ismember(frame, cuts)
        ef = frame;
        nn_blob_size = 0;
        zero_delta_count = 0;
    else  
        if nn_count_frame == 0
            zero_delta_count = zero_delta_count + 1;
            if zero_delta_count >= min_zero_delta_count

                if nn_blob_size > nn_min_blob_size
                    cutsb(frame) = frame + zero_delta_count;
                    ef = frame + zero_delta_count;
                    nn_blob_size = 0;
                    zero_delta_count = 0;
                end
            end
        else
            zero_delta_count = 0;
        end
    end
end

cuts = cuts(~isnan(cuts));
cutsb = cutsb(~isnan(cutsb));

end