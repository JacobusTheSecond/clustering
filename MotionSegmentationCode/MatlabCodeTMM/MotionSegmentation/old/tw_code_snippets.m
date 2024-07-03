% if (options.DEBUG)
%    figure(2);
%    tw_segmentPlot(cuts, subcuts, submots, mot.nframes);
%    print('-depsc', strcat('C:\users\twilli\Desktop\seg_result\', mot.filename, '.eps'));
%    figure(1);
% end

%     if (options.use_mirror_motion)
%         %[ nnidx, nndists ] = tw_makeSymmetric( nnidx, nndists );
%         %[ nnidx_mirror, nndists_mirror ] = tw_makeSymmetric( nnidx_mirror, nndists_mirror );
%         [ nnidx_merge, nndists_merge ] = tw_mergeNN( nnidx, nndists, nnidx_mirror, nndists_mirror, false );     
%         %[ nnidx_merge, nndists_merge ] = tw_makeSymmetric( nnidx_merge, nndists_merge ); 
%     end

%options.frame_offsets = [-10 -5 0 5 10];
%options.frame_offsets = [-9 -6 -3 0 3 6 9];
%options.frame_offsets = [-12 -9 -6 -3 0 3 6 9 12];
%options.frame_offsets = [-8 -6 -4 -2 0 2 4 6 8];
%options.frame_offsets = [-6 -4 -2 0 2 4 6];
%options.frame_offsets = [-5 -3 -1 0 1 3 5];
%options.frame_offsets = [-4 -2 0 2 4];
%options.frame_offsets = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
%options.frame_offsets = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10];
%options.frame_offsets = [-4 -3 -2 -1 0 1 2 3 4];
%options.frame_offsets = [-3 -2 -1 0 1 2 3];
%options.frame_offsets = [-2 -1 0 1 2];
%options.frame_offsets = [-5 -3 0 3 5];
%options.frame_offsets = [-8 -5 0 5 8];

%options.frame_offsets = [-10 -5 -2 0 2 5 10];
%options.frame_offsets = [-10 0 10];
%options.frame_offsets = [-5 0 5];
%options.frame_offsets = [-8 -4 0 4 8];
%options.frame_offsets = [-12 -6 0 6 12]


% check symmetric
% [ nnidx_sym, nndists_sym ] = tw_makeSymmetric( nnidx, nndists );
% a = nnidx-nnidx_sym;
% asum = sum(a(~isnan(a)))
% [ nnidx_mirror_sym, nndists_mirror_sym ] = tw_makeSymmetric( nnidx_mirror, nndists_mirror );
% b = nnidx_mirror-nnidx_mirror_sym;
% bsum = sum(b(~isnan(b)))

%TODO remove
    %[nnidx_filtered ~] = tw_filterRadius(options.radius, nnidx_final, nndists_final);
%     [ y, z, cf, cb ] = tw_countNNperSquare3( nnidx, nndists, options.radius, mot.samplingRate );
%     x = 1:frame_count;
%     x = repmat(x, 5, 1);
%     figure(3);
%     plot(x',y');
%     ylim = get(gca,'ylim');
%     for i=1:numel(cf)
%         line([cf(i) cf(i)], ylim,'color','green','linewidth',2);
%     end
%     figure(4);
%     plot(x',z');
%     size(cb)
%     ylim = get(gca,'ylim');
%     for i=1:numel(cb)
%         line([cb(i) cb(i)], ylim,'color','green','linewidth',2);
%     end
%     figure(5);
%     plot(x(1, :), sum(y));
%     figure(6);
%     plot(x(1, :), sum(z));
%     figure(7);
%     plot(x(1, :), sum(y) + sum(z));
%     figure(1);

%TODO remove
%         [ fmat_foot_l, fmat_foot_r, fmat_head, fmat_hand_l, fmat_hand_r ] = tw_featureSeparation( fmat, 1, frame_count );
%         [ fmat_mirror_foot_l, fmat_mirror_foot_r, fmat_mirror_head, fmat_mirror_hand_l, fmat_mirror_hand_r ] = tw_featureSeparation( fmat_mirror, 1, frame_count );
%         [ nnidx_foot_l, nndists_foot_l, nnidx_mirror_foot_l, nndists_mirror_foot_l] = tw_getMirrorNN( fmat_foot_l, fmat_mirror_foot_l, options.k, options.radius );
%         [ nnidx_merge_foot_l, nndists_merge_foot_l ] = tw_mergeNN(nnidx_foot_l, nndists_foot_l, nnidx_mirror_foot_l, nndists_mirror_foot_l, false);
%         [ nnidx_foot_r, nndists_foot_r, nnidx_mirror_foot_r, nndists_mirror_foot_r ] = tw_getMirrorNN( fmat_foot_r, fmat_mirror_foot_r, options.k, options.radius );
%         [ nnidx_merge_foot_r, nndists_merge_foot_r ] = tw_mergeNN(nnidx_foot_r, nndists_foot_r, nnidx_mirror_foot_r, nndists_mirror_foot_r, false);
%         [ nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head ] = tw_getMirrorNN( fmat_head, fmat_mirror_head, options.k, options.radius );
%         [ nnidx_merge_head, nndists_merge_head ] = tw_mergeNN(nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, false);
%         [ nnidx_hand_l, nndists_hand_l, nnidx_mirror_hand_l, nndists_mirror_hand_l ] = tw_getMirrorNN( fmat_hand_l, fmat_mirror_hand_l, options.k, options.radius );
%         [ nnidx_merge_hand_l, nndists_merge_hand_l ] = tw_mergeNN( nnidx_hand_l, nndists_hand_l, nnidx_mirror_hand_l, nndists_mirror_hand_l, false);
%         [ nnidx_hand_r, nndists_hand_r, nnidx_mirror_hand_r, nndists_mirror_hand_r ] = tw_getMirrorNN( fmat_hand_r, fmat_mirror_hand_r, options.k, options.radius );
%         [ nnidx_merge_hand_r, nndists_merge_hand_r ] = tw_mergeNN(nnidx_hand_r, nndists_hand_r, nnidx_mirror_hand_r, nndists_mirror_hand_r, false);
%         fig = 4;
%         figure(fig);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot_l, nndists_foot_l, nnidx_mirror_foot_l, nndists_mirror_foot_l, nnidx_merge_foot_l, nndists_merge_foot_l); 
%         title('foot_l');
%         figure(fig+1);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot_r, nndists_foot_r, nnidx_mirror_foot_r, nndists_mirror_foot_r, nnidx_merge_foot_r, nndists_merge_foot_r);
%         title('foot_r');
%         figure(fig+2);
%         tw_mergePlot(frame_count, options.radius, nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head); 
%         title('head');
%         figure(fig+3);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand_l, nndists_hand_l, nnidx_mirror_hand_l, nndists_mirror_hand_l, nnidx_merge_hand_l, nndists_merge_hand_l); 
%         title('hand_l');
%         figure(fig+4);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand_r, nndists_hand_r, nnidx_mirror_hand_r, nndists_mirror_hand_r, nnidx_merge_hand_r, nndists_merge_hand_r); 
%         title('hand_r');
%         figure(1);
% 
% 
% %TODO remove
%         [ fmat_foot, fmat_head, fmat_hand ] = tw_featureSeparation3Parts( fmat, 1, frame_count );
%         [ fmat_mirror_foot, fmat_mirror_head, fmat_mirror_hand ] = tw_featureSeparation3Parts( fmat_mirror, 1, frame_count );
%         [ nnidx_foot, nndists_foot, nnidx_mirror_foot, nndists_mirror_foot ] = tw_getMirrorNN( fmat_foot, fmat_mirror_foot, options.k, options.radius );
%         [ nnidx_merge_foot, nndists_merge_foot ] = tw_mergeNN(nnidx_foot, nndists_foot, nnidx_mirror_foot, nndists_mirror_foot, false);
%         [ nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head ] = tw_getMirrorNN( fmat_head, fmat_mirror_head, options.k, options.radius );
%         [ nnidx_merge_head, nndists_merge_head ] = tw_mergeNN(nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, false);
%         [ nnidx_hand, nndists_hand, nnidx_mirror_hand, nndists_mirror_hand ] = tw_getMirrorNN( fmat_hand, fmat_mirror_hand, options.k, options.radius );
%         [ nnidx_merge_hand, nndists_merge_hand ] = tw_mergeNN(nnidx_hand, nndists_hand, nnidx_mirror_hand, nndists_mirror_hand, false);
%         fig = 10;
%         figure(fig);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot, nndists_foot, nnidx_mirror_foot, nndists_mirror_foot, nnidx_merge_foot, nndists_merge_foot);
%         title('foot');
%         figure(fig+1);
%         tw_mergePlot(frame_count, options.radius, nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head); 
%         title('head');
%         figure(fig+2);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand, nndists_hand, nnidx_mirror_hand, nndists_mirror_hand, nnidx_merge_hand, nndists_merge_hand); 
%         title('hand');
%         figure(1);


%         [ fmat_foot_l, fmat_foot_r, fmat_head, fmat_hand_l, fmat_hand_r ] = tw_featureSeparation( fmat, 1, frame_count );
%         [ fmat_mirror_foot_l, fmat_mirror_foot_r, fmat_mirror_head, fmat_mirror_hand_l, fmat_mirror_hand_r ] = tw_featureSeparation( fmat_mirror, 1, frame_count );
%         [ nnidx_foot_l, nndists_foot_l, nnidx_mirror_foot_l, nndists_mirror_foot_l, nnidx_merge_foot_l, nndists_merge_foot_l ] = tw_getAllNN( fmat_foot_l, fmat_mirror_foot_l, options.k, options.radius );
%         [ nnidx_foot_r, nndists_foot_r, nnidx_mirror_foot_r, nndists_mirror_foot_r, nnidx_merge_foot_r, nndists_merge_foot_r ] = tw_getAllNN( fmat_foot_r, fmat_mirror_foot_r, options.k, options.radius );
%         [ nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head ] = tw_getAllNN( fmat_head, fmat_mirror_head, options.k, options.radius );
%         [ nnidx_hand_l, nndists_hand_l, nnidx_mirror_hand_l, nndists_mirror_hand_l, nnidx_merge_hand_l, nndists_merge_hand_l ] = tw_getAllNN( fmat_hand_l, fmat_mirror_hand_l, options.k, options.radius );
%         [ nnidx_hand_r, nndists_hand_r, nnidx_mirror_hand_r, nndists_mirror_hand_r, nnidx_merge_hand_r, nndists_merge_hand_r ] = tw_getAllNN( fmat_hand_r, fmat_mirror_hand_r, options.k, options.radius );
%         fig = 4;
%         figure(fig);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot_l, nndists_foot_l, nnidx_mirror_foot_l, nndists_mirror_foot_l, nnidx_merge_foot_l, nndists_merge_foot_l); 
%         figure(fig+1);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot_r, nndists_foot_r, nnidx_mirror_foot_r, nndists_mirror_foot_r, nnidx_merge_foot_r, nndists_merge_foot_r);
%         figure(fig+2);
%         tw_mergePlot(frame_count, options.radius, nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head); 
%         figure(fig+3);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand_l, nndists_hand_l, nnidx_mirror_hand_l, nndists_mirror_hand_l, nnidx_merge_hand_l, nndists_merge_hand_l); 
%         figure(fig+4);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand_r, nndists_hand_r, nnidx_mirror_hand_r, nndists_mirror_hand_r, nnidx_merge_hand_r, nndists_merge_hand_r); 
%         figure(1);



%         [ fmat_foot, fmat_head, fmat_hand ] = tw_featureSeparation3Parts( fmat, 1, frame_count );
%         [ fmat_mirror_foot, fmat_mirror_head, fmat_mirror_hand ] = tw_featureSeparation3Parts( fmat_mirror, 1, frame_count );
%         [ nnidx_foot, nndists_foot, nnidx_mirror_foot, nndists_mirror_foot, nnidx_merge_foot, nndists_merge_foot ] = tw_getAllNN( fmat_foot, fmat_mirror_foot, options.k, options.radius );
%         [ nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head ] = tw_getAllNN( fmat_head, fmat_mirror_head, options.k, options.radius );
%         [ nnidx_hand, nndists_hand, nnidx_mirror_hand, nndists_mirror_hand, nnidx_merge_hand, nndists_merge_hand ] = tw_getAllNN( fmat_hand, fmat_mirror_hand, options.k, options.radius );
%         fig = 10;
%         figure(fig);
%         tw_mergePlot(frame_count, options.radius, nnidx_foot, nndists_foot, nnidx_mirror_foot, nndists_mirror_foot, nnidx_merge_foot, nndists_merge_foot);
%         figure(fig+1);
%         tw_mergePlot(frame_count, options.radius, nnidx_head, nndists_head, nnidx_mirror_head, nndists_mirror_head, nnidx_merge_head, nndists_merge_head); 
%         figure(fig+2);
%         tw_mergePlot(frame_count, options.radius, nnidx_hand, nndists_hand, nnidx_mirror_hand, nndists_mirror_hand, nnidx_merge_hand, nndists_merge_hand); 
%         figure(1);



%     figure(2);
%     subplot(1,2,1);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_final(1:options.k_s1, :), nndists_final(1:options.k_s1, :));
%     subplot(1,2,2);
%     [nnidx_filtered nndists_filtered] = tw_filterRadius(options.radius * 0.66, nnidx_final, nndists_final);
%     [nnidx_filtered nndists_filtered] = tw_filterDiagonal(15, nnidx_filtered, nndists_filtered);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
%     figure(1);



%     k_filter = floor(size(nnidx_final, 1) / 2);
%     nnidx_k_filter = nnidx_final(1:k_filter, :);
%     nndists_k_filter = nndists_final(1:k_filter, :);
%     p_filter = 0.5;
%     [nnidx_p_filter, nndists_p_filter] = tw_filterNN(nnidx_final, nndists_final, p_filter);
%     figure(2);
%     subplot(1,2,1);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_k_filter, nndists_k_filter);
%     subplot(1,2,1);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_p_filter, nndists_p_filter);
%     figure(1);



%     figure(4);
%     tw_velocityPlot(fmat);
%     figure(1);
%     figure(5);
%     tw_velocityFilterPlot(fmat);
%     figure(1);



%     tw_plotPaths(options, mot, cuts, nnidx_final, nndists_final);
%     [ nnidx_filtered, nndists_filtered ] = tw_possibleStepsFilter( cuts, nnidx_final, nndists_final);
%     figure(2);
%     subplot(1,2,1);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_final, nndists_final);
%     subplot(1,2,2);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_filtered, nndists_filtered);
%     figure(1);



%     figure(15);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx, nndists);
%     %1:options.k_s3
%     %tw_selfSimilarityPlot(frame_count, options.radius, nnidx(1:options.k_s3, :), nndists(1:options.k_s3, :));
%     tw_plotCuts(cuts);
%     tw_plotSubCuts(allsubcuts);
%     figure(1);



%     figure(16);
%     [ nnidx_sym, nndists_sym ] = tw_makeSymmetric( nnidx, nndists );
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_sym, nndists_sym);
%     %1:options.k_s3
%     %tw_selfSimilarityPlot(frame_count, options.radius, nnidx_mirror(1:options.k_s3, :), nndists_mirror(1:options.k_s3, :));
%     tw_plotCuts(cuts);
%     tw_plotSubCuts(allsubcuts);
%     figure(1);
   


%     figure(17);
%     tw_selfSimilarityPlot(frame_count, options.radius, nnidx_mirror, nndists_mirror);
%     %1:options.k_s3
%     %tw_selfSimilarityPlot(frame_count, options.radius, nnidx_mirror(1:options.k_s3, :), nndists_mirror(1:options.k_s3, :));
%     tw_plotCuts(cuts);
%     tw_plotSubCuts(allsubcuts);
%     figure(1);