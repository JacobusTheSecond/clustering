function [cuts,cutsb] = tw_countNNperSquare(nnidx, nn_min_blob_size, min_zero_delta_count, diag_width)
%TW_COUNTNNPERSQUARE Summary of this function goes here
%   Detailed explanation goes here

frame_count = size(nnidx, 2);

cuts = nan(1, frame_count);
cutsb = nan(1, frame_count);

% current number of nearest neighbors in a blob
nn_blob_size = 0;
% minimum number of nearest neighbors in a blob
%nn_min_blob_size = sampling_rate * 4;%128;

% start frame
sf = 1;
% number of zero NN changes
zero_delta_count = 0;
%min_zero_delta_count = floor(sampling_rate / 3);

%% remove main diagonal
%diag_width = floor(sampling_rate);%floor(sampling_rate / 2);%15;
nnidx = tw_filterDiagonal(diag_width, nnidx);

%% find cuts forward 
for i = 1:frame_count
    curidx = nnidx(:, i);

    % remove neighbors that are smaller than sf and greater than i
    curidx(curidx < sf) = nan; %sf
    curidx(curidx > i) = nan;

    nn_i = sum(~isnan(curidx));
    % sum up all nearest neighbors from sf to current frame
    nn_blob_size = nn_blob_size + nn_i;
    
    if nn_i == 0
        zero_delta_count = zero_delta_count + 1;
        
        % if the nearest neighbor count doesn't change 
        % max_zero_delta_count times in a row
        if zero_delta_count >= min_zero_delta_count
            
            % If the blob size is large enough, then make a cut.
            % If any frame has more than minblobsize NN
            if nn_blob_size > nn_min_blob_size
                cuts(i) = i - zero_delta_count;
                sf = i - zero_delta_count;
                nn_blob_size = 0;
                zero_delta_count = 0;
            end
        end
    else
        zero_delta_count = 0;
    end
end

%% find cuts backwards

% end frame
ef = frame_count;
nn_blob_size = 0;
zero_delta_count=0;

for i = frame_count:-1:1
    curidx = nnidx(:, i);

    % remove neighbors that are smaller than i and greater than ef
    curidx(curidx < i) = nan;
    curidx(curidx > ef) = nan; %ef
    
    nn_i = sum(~isnan(curidx));
    % sum up all nearest neighbors from ef to current frame
    nn_blob_size = nn_blob_size + nn_i;
 
    % if frame is a cut then ignore the current blob
    if ismember(i, cuts)
        ef = i;
        nn_blob_size = 0;
        zero_delta_count = 0;
    else  
        if nn_i == 0
            zero_delta_count = zero_delta_count + 1;
            if zero_delta_count >= min_zero_delta_count

                if nn_blob_size > nn_min_blob_size
                    cutsb(i) = i + zero_delta_count;
                    ef = i + zero_delta_count;
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