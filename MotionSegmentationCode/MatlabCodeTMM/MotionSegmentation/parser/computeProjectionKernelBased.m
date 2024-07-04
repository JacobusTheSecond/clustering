function [pdata,meta]=computeProjectionKernelBased(data,varargin)

options.DEBUG         = false;
options.VIZPROJECTION = false;
options.pcaThreshold  = 0.98;
options.k             = 32;
options.kernel        = 'gaussian';
% options.kernel       = 'sqrt';

switch nargin
    case 2
        options = mergeOptions(varargin{1},options);
end

tic;
t1 = toc;

% reduce to pca space
dataorg = data;

[coeff,score,latent] = princomp2(data');
ncomps = find(cumsum(latent)/sum(latent)>=options.pcaThreshold,1,'first');

ncomps = max(ncomps,3);

meanVec = mean(data,2);
data    = score(:,1:ncomps)';

nframes = size(data,2);
ndims   = size(data,1);


% get knn per frame:
act_k         = max(options.k+1,ncomps+1); % make sure we have enough variables for the optimization
options.act_k = min(nframes-1,act_k);

fprintf('Performing %i knn searches (%iNN) in %i dimensions... ',nframes,options.act_k-1,ndims);

t2 = toc;

handle = ann_buildTree(data);
[oridx,ordists] = ann_queryTree(handle,data,options.act_k);
ann_cleanup();

oridx   = oridx(2:end,:);
ordists = ordists(2:end,:);

fprintf('done!\n')

t3 = toc;

% get local directions:
fprintf('compute directions... ')
dirs = computeDirections_local(data);
fprintf('done!\n')
% 
options_optim = optimset( 'Display','none', ...
                          'MaxFunEvals',2^13, ...
                          'TolX',0.00001, ...
                          'TolFun',0.00001, ...
                          'Algorithm','levenberg-marquardt');

% copy data for later results
pdata = zeros(size(data));

fprintf('Moving points:             ')
parfor f = 1:nframes

    fprintf('\b\b\b\b\b\b\b\b\b\b\b')
    fprintf('%05i/%05i',f,nframes);
    
    dir = dirs(:,f);
    
    R           = getRotationMat_local(dir,ndims);
    knnweights	= costs2weights_local(ordists(:,f));
    
% % %     if options.DEBUG
% % %         figure();
% % %         hold on;
% % %         plot3(data(1,f),data(2,f),data(3,f),'*','color','red');
% % %         plot3(data(1,oridx(:,f)),data(2,oridx(:,f)),data(3,oridx(:,f)),'.');
% % %         line([data(1,f) data(1,f)+R(1,1)],[data(2,f) data(2,f)+R(2,1)],[data(3,f) data(3,f)+R(3,1)],'color','red');
% % %         line([data(1,f) data(1,f)+R(1,2)],[data(2,f) data(2,f)+R(2,2)],[data(3,f) data(3,f)+R(3,2)],'color','green');
% % %         line([data(1,f) data(1,f)+R(1,3)],[data(2,f) data(2,f)+R(2,3)],[data(3,f) data(3,f)+R(3,3)],'color','blue');
% % %          line([data(1,f) data(1,f)+dir(1)],[data(2,f) data(2,f)+dir(2)],[data(3,f) data(3,f)+dir(3)],'color','red','linewidth',3);
% % %         grid on;
% % %         axis equal;
% % %     end
    
    if options.DEBUG && mod(f,16)==0;
       
        figure();
        set(gcf,'position',[70 70 1000 1000]);
        
        xmin = min(data(1,oridx(:,f)));
        xmax = max(data(1,oridx(:,f)));
        ymin = min(data(2,oridx(:,f)));
        ymax = max(data(2,oridx(:,f)));
        
        xs = xmin:0.02:xmax;
        ys = ymin:0.02:ymax;
        
        r_g = zeros(numel(xs),numel(ys));
        r_s = zeros(numel(xs),numel(ys));
        
        for i=1:numel(xs)
            for j=1:numel(ys)
                
                dists = sqrt(sum((repmat([xs(i);ys(j)],1,options.act_k-1) - data(1:2,oridx(:,f)) ).^2))';
                r_g(i,j) = computeGaussian_local(dists, ones(size(oridx,1),1));                
%                 r_s(i,j) = computePrior_local(dists, ones(size(oridx,1),1),0.5);
            end
        end

%         subplot(1,2,1);
        imagesc(xs,ys,r_g');
        cmap = flipud(hot(128));
        colormap(cmap);
        
        hold on;
        
        curidx = sort([oridx(:,f);f]);
        idxD = diff(curidx);
        idxJ = [1; find(idxD>1)+1; numel(curidx)];
        
        for seq = 1:numel(idxJ)-1
            plot(data(1,curidx(idxJ(seq):idxJ(seq+1)-1)), ...
                 data(2,curidx(idxJ(seq):idxJ(seq+1)-1)), ...
                 '.-','linewidth',1,'color','green');
        end
        plot([data(1,f)-dir(1)*10 data(1,f)+dir(1)*10], ...
             [data(2,f)-dir(2)*10 data(2,f)+dir(2)*10], ...
              '-','linewidth',3,'color','blue')
        
        plot(data(1,f),data(2,f),'o','linewidth',3,'color','red');
        colorbar;
        grid on;
        title(['Gaussian density estimation Frame ' num2str(f)] );
%         subplot(1,2,2);
%         imagesc(xs,ys,r_s');
%         colormap hot;
%         hold on;
%         plot(data(1,oridx(:,f)),data(2,oridx(:,f)),'+','linewidth',3,'color','green');
%         colorbar;
%         grid on;
%         title('sqrt-kernel density estimation')
        
    end

    knndata	= data(:,oridx(:,f));
    odata	 = data(:,f);
    
    % start value (original point)
    V0 = zeros(ndims-1,1);

    V = lsqnonlin(   @(V)objfun_local(  V,odata,knndata,knnweights,R,options), ...
                     V0,[],[],options_optim );
                                
    % move point to final position:
    V = [0;V];
    pdata(:,f) = R * V + data(:,f);
%     pdata(:,f) = V + data(:,f);
    
% %     if options.DEBUG && mod(f,16)==0;
% %         plot(pdata(1,f),pdata(2,f),'+','linewidth',3,'color','red');
% %     end

% % %     if options.DEBUG
% % %         plot3(pdata(1,f),pdata(2,f),pdata(3,f),'*','color','red','linewidth',3);
% % %     end
    
end

fprintf('done!\n\n');

%pca back projection:
pdata = coeff(:,1:ncomps) * pdata + meanVec(:,ones(1,nframes)); 

t4 = toc;

if options.VIZPROJECTION; %options.DEBUG 
    
%     close all;
    
% % %     figure()
% % %     subplot(3,1,1);
% % %     plot(dataorg(1:3,:))
% % %     plot(dataorg(1:3,:)')
% % %     grid on
% % %     subplot(3,1,2);
% % %     plot(pdata(1:3,:)')
% % %     grid on
% % %     subplot(3,1,3);
% % %     plot(data(1:2,:)')
% % %     grid on
% % %     
% % % % %     figure()
% % % % %     subplot(1,2,1);
% % % % %     plot3(dataorg(1,:),dataorg(2,:),dataorg(3,:),'.-')
% % % % %     grid on; axis equal;
% % % % %     subplot(1,2,2);
% % % % %     plot3(pdata(1,:),pdata(2,:),pdata(3,:),'.-')
% % % % %     grid on; axis equal;
% % %     

%     pdata = [pdata(1:end/2,:) pdata(end/2+1:end,:)];
%     dataorg = [dataorg(1:end/2,:) dataorg(end/2+1:end,:)];
% 
%     myf = figure();
%     set(myf,'position',[50 50 1500 800]);
%     subplot(1,2,1);
%     [coeff,score,latent] = princomp2(dataorg');
%     plot3(score(1:end/2,1),score(1:end/2,2),score(1:end/2,3),'.-','color','red');
%     hold on
%     plot3(score(end/2:end,1),score(end/2:end,2),score(end/2:end,3),'.-','color','blue');
%     grid on; axis equal;
%     
%     for i=1:1:size(score,1)/2
%         line([score(i,1) score(i+size(score,1)/2,1)],[score(i,2) score(i+size(score,1)/2,2)],[score(i,3) score(i+size(score,1)/2,3)],'color',[0.5 0.5 0.5])
%     end
%     
%     
%     subplot(1,2,2);
%     [coeff,score,latent] = princomp2(pdata');
%     plot3(score(1:end/2,1),score(1:end/2,2),score(1:end/2,3),'.-','color','red');
%     hold on
%     plot3(score(end/2:end,1),score(end/2:end,2),score(end/2:end,3),'.-','color','blue');
%     grid on; axis equal;
%     
%     for i=1:1:size(score,1)/2
%         line([score(i,1) score(i+size(score,1)/2,1)],[score(i,2) score(i+size(score,1)/2,2)],[score(i,3) score(i+size(score,1)/2,3)],'color',[0.5 0.5 0.5])
%     end
%     
%     
%     drawnow();
    

    figure()
    dataviz = [dataorg pdata];
    [~,score,~] = princomp2(dataviz');
    subplot(1,2,1);
    plot3(score(1:end/2,1),score(1:end/2,2),score(1:end/2,3),'.-','color','blue');
    grid on; axis equal;
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    zlims = get(gca,'zlim');
    
% % %     hold on
% % %     plot3(score([395 405 415],1),score([395 405 415],2),score([395 405 415],3),'o','color','red','linewidth',3);
% % %     plot3(score([495 505 515],1),score([495 505 515],2),score([495 505 515],3),'o','color','green','linewidth',3);
% % %     plot3(score([750 760 770],1),score([750 760 770],2),score([750 760 770],3),'o','color','cyan','linewidth',3);
% % %     hold off
    
    s2=subplot(1,2,2);
    plot3(score(end/2+1:end,1),score(end/2+1:end,2),score(end/2+1:end,3),'.-','color','blue');
    grid on; axis equal;
    set(s2,'xlim',xlims);
    set(s2,'ylim',ylims);
    set(s2,'zlim',zlims);
    
% % %     hold on
% % %     plot3(score(end/2+[395 405 415],1),score(end/2+[395 405 415],2),score(end/2+[395 405 415],3),'o','color','red','linewidth',3);
% % %     plot3(score(end/2+[495 505 515],1),score(end/2+[495 505 515],2),score(end/2+[495 505 515],3),'o','color','green','linewidth',3);
% % %     plot3(score(end/2+[750 760 770],1),score(end/2+[750 760 770],2),score(end/2+[750 760 770],3),'o','color','cyan','linewidth',3);
% % %     hold off
    
    set(gcf,'position',[70 70 1800 1000]);
    

    
    

% %     myf = figure();
% %     set(myf,'position',[50 50 1500 800]);
% %     subplot(1,2,1);
% %     [coeff,score,latent] = princomp2(dataorg');
% %     plot3(score(:,1),score(:,2),score(:,3),'.-','color','blue');
% %     grid on; axis equal;
% %   
% %     
% %     subplot(1,2,2);
% %     [coeff,score,latent] = princomp2(pdata');
% %     plot3(score(:,1),score(:,2),score(:,3),'.-','color','blue');
% %   
% %     grid on; axis equal;
 
    
    drawnow();



end

meta.time     = t2-t1 + t4-t3;
meta.dims     = ncomps;
meta.accuracy = options.pcaThreshold;
meta.knn      = options.k;
meta.kernel   = options.kernel;

end

function F = objfun_local(V,odata,knndata,knnweights,R,options)

V = [0;V];
ndata = R * V + odata;
% ndata = V + odata;

knndists = sqrt(sum((knndata-ndata(:,ones(1,size(knndata,2)))).^2))';
% knndists = pdist2(knndata',ndata');

switch options.kernel
    case 'sqrt'
        F = computePrior_local( ...
                                knndists, knnweights, 0.5);
    case 'gaussian'
        F = computeGaussian_local( ...
                                knndists, knnweights);
    otherwise
        error('unknown kernel!')
end
                    
% % % if options.DEBUG
% % %     plot3(ndata(1),ndata(2),ndata(3),'.','color','green');
% % % end

end

function Q = getRotationMat_local(dir,dim)

    A = rand(dim);
    A(:,1) = dir;
    
    A = normalizeColumns(A);
    
    [Q,~]=qr(A);

end

function dirs = computeDirections_local(data)

filterSize = 5;

nframes = size(data,2);

jointVelocities     = zeros(size(data));

if size(data,2)>1
    % padding
    data = [	3*data(:,1)-2*data(:,2),...
                2*data(:,1)-data(:,2),...
                  data,...
                2*data(:,end)-data(:,end-1),...
                3*data(:,end)-2*data(:,end-1)];

    % 5-point derivation
    weights = [1 -8 0 8 -1];
    for frame = 3:nframes+2
        jointVelocities(:,frame-2) = ...
                  weights(1) * data(:,frame-2) ...
                + weights(2) * data(:,frame-1) ...
                + weights(3) * data(:,frame)   ...
                + weights(4) * data(:,frame+1) ...
                + weights(5) * data(:,frame+2);
    end  
    jointVelocities = jointVelocities / (12);

    % filtering velocities
    jointVelocities = filterTimeline(jointVelocities,filterSize,'bin');
    
end

dirs = normalizeColumns(jointVelocities);

end

function e = computePrior_local(diff,poseWeights,exponent)
    e = diff.^exponent; 
    e = sum(e.* poseWeights);
end

function e = computeGaussian_local(diff,poseWeights)

    sigma  =  1/(mean(abs(diff)));
    factor = (1/(sigma*sqrt(2*pi)));
    
    e = factor * exp(-(1/(2*sigma^2)) * diff.^2);

%     sigma = std(diff);
%     mu    = mean(diff);
% 
%     e = 1/(sigma*sqrt(2*pi))*exp(-(diff-mu).^2/(2*sigma*sigma));
    
    e = factor * size(diff,1) - sum(e.* poseWeights);
%      e = sum(e.* poseWeights);
end

function weights = costs2weights_local(costs)
    if ~any(costs~=0)
        weights = ones(size(costs))/numel(costs);
    else
        weights = costs-min(costs);
        weights = weights/max(weights);
        weights = 1-weights;
    end
end
