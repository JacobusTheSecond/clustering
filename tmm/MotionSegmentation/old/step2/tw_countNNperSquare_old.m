function [ nns,cuts,cutsb ] = tw_countNNperSquare_old( nnidx )
%TW_COUNTNNPERSQUARE_OLD Summary of this function goes here
%   Detailed explanation goes here

DEBUG = false;

%facial params
% % % mdremoval = 'search';
% % % minblobsize = 512;%256;

% % % mocap params
mdremoval = 'fixed';
minblobsize = 128;

nns  = nan(1,size(nnidx,2));

% remove main diag!
tic;
switch mdremoval
    case 'fixed'
        diagwidth = 15;
        for i=1:size(nnidx,2)

            nnidx(abs(nnidx(:,i)-i)<diagwidth,i)=nan;

        end
    case 'search'
        for i=1:size(nnidx,2)
            curidx = nnidx(:,i);
            
            [curidx] = sort(curidx);
            
            idxdiff = diff(curidx);
            ind = find(curidx==i);
            
            blocks = find(idxdiff>1);
            sid = find(blocks<=ind,1,'last');
            eid = find(blocks>=ind,1,'first');
            
            if isempty(sid)
                sidx = 1;
            else
                sidx = blocks(sid);
            end
            
            if isempty(eid)
                eidx = size(curidx,1);
            else
                eidx = blocks(eid);
            end
            
            curidx(sidx:eidx) = nan;

            nnidx(:,i)=curidx;

        end
    otherwise
        error('Unknown method to remove main diag!');
end

% % fprintf('Removing main diag took: %2.2f sec.\n',toc);

if DEBUG
    
    immat = zeros(size(nnidx,2));
    for i=1:size(nnidx,2)
        immat(i,nnidx(~isnan(nnidx(:,i)),i))=1;
    end

    figure();
    colormap gray;
    matfigure = imagescnan([],[],immat',isnan(immat'),'b');
    drawnow()
    
end



sf      = 1;
curZ    = 0;
lastcut = 1;

% % % for c=1:36
% % %     fprintf('.');
% % % end

if DEBUG 
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    lhvi  = line([1 1], ylim,'color','green');
    lhhi  = line(xlim, [1 1],'color','green');
    lhvsf  = line([1 1], ylim,'color','red');
    lhhsf  = line(xlim, [1 1],'color','red');
end


for i=2:size(nnidx,2)
    curidx    = nnidx(:,sf:i);

    curidx(curidx<sf)= nan;
    curidx(curidx>i) = nan;
    
    nns(i) = sum(sum(~isnan(curidx)));

	if DEBUG
       
        set(lhvi,'xdata',[i i]);
        set(lhhi,'ydata',[i i]);
        set(lhvsf,'xdata',[sf sf]);
        set(lhhsf,'ydata',[sf sf]);
        drawnow()
        title(num2str(nns(i)));
    end
    
    
    % ausgaben:
% % %     for c=1:36
% % %         fprintf('\b');
% % %     end
% % %     fprintf('sf = %04i, ef = %04i, nnz = %08i',sf,i,nns(i));
    
    cd = nns(i) - nns(i-1);
    
    if cd == 0
        curZ = curZ+1;
        if curZ>=10
            
            if any(nns(lastcut:i)>minblobsize)
                lastcut = i-curZ;
                sf = i;
                curZ = 0;
            else
                curZ=0;
            end
        end
    else
        curZ=0;
    end  
   
end

% % % fprintf('\n')
% % % for c=1:36
% % %         fprintf('.');
% % % end
nnsb  = nan(1,size(nnidx,2));
ef    = size(nnidx,2);
lastcut = size(nnidx,2);

if DEBUG 
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    lhvi  = line([1 1], ylim,'color','green');
    lhhi  = line(xlim, [1 1],'color','green');
    lhvef  = line([1 1], ylim,'color','red');
    lhhef  = line(xlim, [1 1],'color','red');
end

[~,cuts]  = findpeaks(nns, 'MINPEAKDISTANCE',30,'MINPEAKHEIGHT',minblobsize);
% cuts=[];
curZ=0;
for i=size(nnidx,2)-1:-1:1
    curidx    = nnidx(:,i:ef);

    curidx(curidx<i)  = nan;
    curidx(curidx>ef) = nan;
    
    nnsb(i) = sum(sum(~isnan(curidx)));

    if DEBUG
       
        set(lhvi,'xdata',[i i]);
        set(lhhi,'ydata',[i i]);
        set(lhvef,'xdata',[ef ef]);
        set(lhhef,'ydata',[ef ef]);
        drawnow()
        title(num2str(nnsb(i)));
    end
    
    
    % ausgaben:
% % %     for c=1:36
% % %         fprintf('\b');
% % %     end
% % %     fprintf('sf = %04i, ef = %04i, nnz = %08i',i,ef,nnsb(i));
    
    cd = nnsb(i) - nnsb(i+1);
    
    if ismember(i,cuts)
        lastcut = i;
        ef = i;
    else
    
        if cd == 0
            curZ = curZ+1;
            if curZ>=10

                if any(nnsb(i:lastcut)>minblobsize)

                    curZ = 0;

                    lastcut = i+curZ;
                    ef = i+curZ;
                else
                    curZ=0;
                end
            end
        else
            curZ=0;
        end
    end
end



[~,cutsb] = findpeaks(nnsb,'MINPEAKDISTANCE',30,'MINPEAKHEIGHT',minblobsize);

% add curZ width:
cutsb = cutsb+10;

% % % fprintf('\n');

if DEBUG
    figure()
    subplot(2,1,1)
    plot(nns)
    hold all
    grid on
    ylim = get(gca,'ylim');

    cuts = cuts+1;
    for i=1:numel(cuts)
        line([cuts(i) cuts(i)], ylim,'color','green');
    end

    subplot(2,1,2)
    plot(nnsb)
    grid on
    ylim = get(gca,'ylim');

    cutsb = cutsb-1;
    for i=1:numel(cutsb)
        line([cutsb(i) cutsb(i)], ylim,'color','red');
    end
    
%     xlim = get(matfigure,'xlim');
%     ylim = get(matfigure,'ylim');
%     for i=1:numel(cuts)
%         line([cuts(i) cuts(i)], ylim,'color','green');
%         line(xlim, [cuts(i) cuts(i)],'color','green');
%     end
% 
%     for i=1:numel(cutsb)
%         line([cutsb(i) cutsb(i)], ylim,'color','red');
%         line(xlim, [cutsb(i) cutsb(i)],'color','red');
%     end
    
end

end
