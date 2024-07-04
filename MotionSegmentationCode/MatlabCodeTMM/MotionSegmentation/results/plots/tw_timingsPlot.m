function [timings_table] = tw_timingsPlot( alltimings, frame_counts, names, all_name )
%TW_TIMINGSPLOT Summary of this function goes here
%   Detailed explanation goes here

    alltimings_count = numel(alltimings);
    timings_names = fieldnames(alltimings{1});
    timings_count = numel(timings_names);

    A = zeros(alltimings_count + 1, timings_count);
    
    for row=1:alltimings_count
        timing = alltimings{row};
        for column=1:timings_count
            A(row, column) = timing.(timings_names{column});
        end
    end
    
    all_frames = floor(sum(frame_counts) / numel(frame_counts));
    weights = frame_counts / sum(frame_counts);
    A_weights = repmat(weights(:), 1, timings_count);
    % normal average
    %A(end, :) = (sum(A, 1) / alltimings_count);
    % weighted average
    A(end, :) = sum(A(1:end-1, :) .* A_weights, 1);
    
    % set legend names
    legend_names = cell(size(timings_names));
    for i=1:numel(timings_names)
        legend_names{i} = strrep(timings_names{i}, '_', ' ');
    end
    
    % set trial names
    trial_names = cell(size(names));
    for i=1:numel(trial_names)
        trial_names(i) = strcat(names{i}, {' | '}, num2str(frame_counts(i)));
    end
    trial_names(numel(trial_names) + 1) = strcat(all_name, {' | '}, num2str(all_frames));
    
    A = fliplr(A);
    legend_names = flipud(legend_names);
        
    ax = gca;
    h = bar(A, 'stacked');
    order=numel(legend_names):-1:1;
    legend(h(order), legend_names{order});
    xlim([0.5 (numel(trial_names)+0.5)]);
    set(ax, 'XTick', 1:numel(trial_names));
    set(ax, 'YTick', 0:100);
    set(ax, 'xticklabel', trial_names);
    xlabel('Trial | Frame Count');
    ylabel('Execution Time [s]');
    set(ax, 'XGrid', 'off');
    set(ax, 'YGrid', 'on');
    colormap(jet(256));
    
    hh = rotateXLabels(ax, 45);
    set(hh, 'interpreter', 'none');
    
    % create timings table
    tt_row_count = 1 + numel(trial_names);
    tt_column_count = 1 + numel(legend_names) + 1;
    timings_table = cell(tt_row_count, tt_column_count);
    
    timings_table{1,1} = 'Trial | Frame Count';
    timings_table(2:numel(trial_names)+1, 1) = trial_names;
    timings_table(1, numel(legend_names)+1:-1:2) = cellstr(strcat(legend_names, ' [s]'));
    timings_table{1, numel(legend_names)+2} = 'Entire Trial [s]';
    A_str = cell(size(A));
    for i = 1:size(A, 2)
       A_str(:, i) = cellstr(num2str(A(:, i), '%.2f'));   
    end
    timings_table(2:end, numel(legend_names)+1:-1:2) = A_str;
    
    entire_trial_time = sum(A, 2);
    entire_trial_time = entire_trial_time(:);
    timings_table(2:end, numel(legend_names)+2) = cellstr(num2str(entire_trial_time, '%.2f'));
end

