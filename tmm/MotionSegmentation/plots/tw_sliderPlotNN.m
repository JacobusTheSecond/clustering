function [] = tw_sliderPlotNN(nndists, radius)
%TW_SLIDERPLOT Summary of this function goes here
%   Detailed explanation goes here

%     dbstop in title
%     title('abc')

    k = size(nndists, 1);
    frame_count = size(nndists, 2);

    figure(2);
    set(gcf, 'position',[100 100 600 600]);
    
    frame = 1;
    nndist = nndists(:,frame);
    nndist = nndist(~isinf(nndist));
    plot(1:numel(nndist), nndist, 'b');
    legend(num2str(frame));
    axis manual;
    axis([1 k 0 radius]);
    
    uicontrol('style','slide',...
                     'position',[80 10 470 20],...
                     'min',1,'max',frame_count,'val',1,...
                     'callback',{@sl_call, nndists, radius});
end

function [] = sl_call(varargin)
    % Callback for the slider.
    [h, nndists, radius] = varargin{[1,3,4]};  % calling handle and data structure.
    
    cla
    k = size(nndists, 1);
    
    frame = round(get(h,'value'));    
    nndist = nndists(:,frame);
    nndist = nndist(~isinf(nndist));
    plot(1:numel(nndist), nndist, 'b');
    legend(num2str(frame));
    axis manual;
    axis([1 k 0 radius]);
end