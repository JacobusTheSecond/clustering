function [strict_results, tolerant_results] = tw_evaluation( range, options_list, record_names, base_path, skel_file, file_prefix, file_suffix )
%TW_EVALUATION Summary of this function goes here
%   Detailed explanation goes here

    range_count = numel(range);
    options_list_count = numel(options_list);
    names = cell(1, range_count);
    frame_counts = zeros(range_count, options_list_count);
    strict_results = zeros(range_count, options_list_count);
    tolerant_results = zeros(range_count, options_list_count);
    
    [ gt, annotation ] = tw_gt_86();
    for options_list_idx=1:options_list_count
        options = options_list{options_list_idx};
        
        for i=range
           gt_data = gt{i};
           annotation_data = annotation{i};
           [ gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color ] = tw_gt_to_segments(gt_data, annotation_data);
           [ skel, mot, name ] = tw_get_mot_by_index(i, base_path, skel_file, file_prefix, file_suffix);
           [mot, mots, submots, comps, cuts, subcuts, subcuts_main, subcuts_mirror, timings] = tw_segmentation(skel, mot, options);
           frame_count = mot.nframes;
           [ cut_segment_start, cut_segment_end ] = tw_cuts_to_segments( cuts, frame_count );
           [ strict, tolerant, cut_segment_to_gt_segment, cut_segment_color ] = tw_segment_match( gt_segment_start, gt_segment_end, gt_is_transition, gt_segment_color, cut_segment_start, cut_segment_end);

           names{i} = name;
           frame_counts(i, options_list_idx) = frame_count;
           strict_results(i, options_list_idx) = strict;
           tolerant_results(i, options_list_idx) = tolerant;

           figure(10 + i);
           subplot(2,1,1);
           tw_plot_segmentation( gt_segment_start, gt_segment_end, gt_segment_color )
           subplot(2,1,2);
           tw_plot_segmentation( cut_segment_start, cut_segment_end, cut_segment_color )
        end
    end
    
    [ strict_old, tolerant_old, old_record_names ] = tw_old_evaluation();

    record_names_ext = cell(1, numel(record_names) + numel(old_record_names));
    for j = 1:numel(record_names)
        record_names_ext{j} = record_names{j};
    end
    for j = 1:numel(old_record_names)
        record_names_ext{numel(record_names) + j} = old_record_names{j};
    end
    figure(101);
    strict_old_range = strict_old(range, :);
    strict_records_count = size(strict_results, 2);
    strict_old_records_count = size(strict_old_range, 2);
    strict_ext = zeros(range_count, strict_records_count + strict_old_records_count);
    strict_ext(:, 1:strict_records_count) = strict_results;  
    strict_ext(:, (strict_records_count + 1):end) = strict_old_range;
    tw_accuracyPlot( strict_ext, 'Strict Evaluation', names, record_names_ext, lines(size(strict_ext, 2)) );
    figure(102);
    tolerant_old_range = tolerant_old(range, :);
    tolerant_records_count = size(tolerant_results, 2);
    tolerant_old_records_count = size(tolerant_old_range, 2);
    tolerant_ext = zeros(range_count, tolerant_records_count + tolerant_old_records_count);
    tolerant_ext(:, 1:tolerant_records_count) = tolerant_results;  
    tolerant_ext(:, (tolerant_records_count + 1):end) = tolerant_old_range;
    tw_accuracyPlot( tolerant_ext, 'Tolerant Evaluation', names, record_names_ext, lines(size(tolerant_ext, 2)) );
end
