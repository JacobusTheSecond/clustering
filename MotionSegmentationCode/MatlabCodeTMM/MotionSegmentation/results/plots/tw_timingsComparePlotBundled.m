function [timings_table] = tw_timingsComparePlotBundled( first_alltimings, second_alltimings, frame_counts, first_names, second_names, first_all_name, second_all_name, title_name )
%TW_TIMINGSCOMPAREPLOT Summary of this function goes here
%   Detailed explanation goes here

    alltimings_count = numel(first_alltimings) * 2;
    timings_names = fieldnames(first_alltimings{1});
    timings_count = numel(timings_names);
    trial_names = cell(1, alltimings_count);
    
    A = zeros(alltimings_count + 2, timings_count);
    
    type = true;    
    for row=1:alltimings_count
        if type
            index = floor(row / 2) + 1;
            timing = first_alltimings{index};
            trial_names(row) = strcat(first_names{index}, {' | '}, num2str(frame_counts(index)));
        else
            index = floor(row / 2);
            timing = second_alltimings{index};
            trial_names(row) = strcat(second_names{index}, {' | '}, num2str(frame_counts(index)));
        end
        
        type = ~type;
        
        for column=1:timings_count
            A(row, column) = timing.(timings_names{column});
        end
    end
    
    all_frames = floor(sum(frame_counts) / numel(frame_counts));
    weights = frame_counts / sum(frame_counts);
    A_weights = repmat(weights(:), 1, timings_count);
    % weighted average
    A(end-1, :) = sum(A(1:2:end-2, :) .* A_weights, 1);
    A(end, :) = sum(A(2:2:end-2, :) .* A_weights, 1);
    
    % set trial names
    trial_names(numel(trial_names) + 1) = strcat(first_all_name, {' | '}, num2str(all_frames));
    trial_names(numel(trial_names) + 1) = strcat(second_all_name, {' | '}, num2str(all_frames));
       
    % set legend names
    legend_names = cell(size(timings_names));
    for i=1:numel(timings_names)
        legend_names{i} = strrep(timings_names{i}, '_', ' ');
    end
    
    A = fliplr(A);
    legend_names = flipud(legend_names);
        
    ax = gca;
    h = bar(A, 'stacked');
    order=numel(legend_names):-1:1;
    legend(h(order), legend_names{order});
    xlim([0.5 (numel(trial_names)+0.5)]);
    set(ax, 'XTick', 1:numel(trial_names));
    set(ax, 'YTick', 0:5:180);
    set(ax, 'xticklabel', trial_names);
    xlabel('Trial | Frame Count');
    ylabel('Execution Time [s]');
    set(ax, 'XGrid', 'off');
    set(ax, 'YGrid', 'on');
    colormap(jet(256));
    title(title_name);
    
    % draw lines    
    trial_count = numel(trial_names);
    if (trial_count > 2)
        lim = ylim;
        line_range=2.5:2:trial_count-1.5;
        for i = 1:numel(line_range)
            line([line_range(i) line_range(i)], lim, 'Color', 'k', 'LineWidth', 2);
        end
    end
    
    hh = rotateXLabels(ax, 45);
    set(hh, 'interpreter', 'none');
    
    % create timings table
    tt_row_count = 1 + numel(trial_names);
    tt_column_count = 1 + numel(legend_names) + 1 + 1;
    timings_table = cell(tt_row_count, tt_column_count);
    
    timings_table{1,1} = 'Trial | Frame Count';
    timings_table(2:numel(trial_names)+1, 1) = trial_names;
    timings_table(1, numel(legend_names)+1:-1:2) = cellstr(strcat(legend_names, ' [s]'));
    timings_table{1, numel(legend_names)+2} = 'Entire Trial [s]';
    timings_table{1, numel(legend_names)+3} = 'Speedup per Trial';
    A_str = cell(size(A));
    for i = 1:size(A, 2)
       A_str(:, i) = cellstr(num2str(A(:, i), '%.2f'));   
    end
    timings_table(2:end, numel(legend_names)+1:-1:2) = A_str;
    
    entire_trial_time = sum(A, 2);
    entire_trial_time = entire_trial_time(:);
    timings_table(2:end, numel(legend_names)+2) = cellstr(num2str(entire_trial_time, '%.2f'));
    
    speedup_increase = entire_trial_time(2:2:end) ./ entire_trial_time(1:2:end);
    speedup_increase_strs = num2str(speedup_increase, '%.2f');
    timings_table(2:2:end, numel(legend_names)+3) = cellstr(speedup_increase_strs);
    timings_table(3:2:end, numel(legend_names)+3) = {' '};
end

