function [A_all] = tw_accuracyPlot( A, title_name, names, record_names, colors, frame_counts )
%TW_ACCURACYPLOT Summary of this function goes here
%   Detailed explanation goes here

    allnames = names;
    allnames{numel(names) + 1} = 'Average';
    
    A_all = zeros(size(A, 1) + 1, size(A, 2));
    A_all(1:size(A, 1), :) = A;
    
    A_all(end, :) = A'*frame_counts(:,1)./sum(frame_counts(:,1));
        
    
    %A_all(end, :) = sum(A) / size(A, 1);
    
    ax = gca;
    bar(A_all(:,1:end));
    set(ax, 'xticklabel', allnames);
    xlabel('Trial');
    ylabel('Accuracy');
    set(ax, 'Ylim', [0.4 1]);
    set(ax, 'XGrid', 'off');
    set(ax, 'YGrid', 'on');
    title(title_name);
    legend(record_names(1:end));
    colormap( colors(1:size(A_all, 2), :) );
end

