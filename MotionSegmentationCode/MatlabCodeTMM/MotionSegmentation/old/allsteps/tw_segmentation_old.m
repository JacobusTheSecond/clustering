function [mots,submots,comps, allsubcuts, timings] = tw_segmentation_old(skel, mot, varargin)
%TW_SEGMENTATION_OLD Summary of this function goes here
%   Detailed explanation goes here

%%
%preliminaries
close all;
dopts.DEBUG    = false;
%
% R FOOT 1,2,3
% L FOOT 4,5,6
% HEAD 7,8,9
% R HAND 10,11,12
% L HAND 13,14,15
%
% relativ positions from root
% y == height
%

dopts.feature_set = 'e15_flex';
%dopts.frame_offsets = [-10 -5 0 5 10];
dopts.frame_offsets = [-5 0 5];
dopts.k = 600;
dopts.generalized_radius = 30.1;
dopts.generalized_radius_activities = 27;

dopts.minrad    = 16; 
dopts.maxrad    = 128;

%  Parameters for segmentation of the facial data from Jun
% dopts.DEBUG    = true;
% dopts.datatype = 'facial';
% dopts.fset     = 'e15_45';
% dopts.k        = floor(mot.nframes/16);
% dopts.radius   = 0.1;
% dopts.radius_s1= 1;
% dopts.k_s1     = floor(dopts.k/2);
% dopts.radius_s3= 1;
% dopts.k_s3     = floor(dopts.k/2);
% dopts.minrad   = 4; 
% dopts.maxrad   = 64;


dopts.minwidth = 0.25*mot.samplingRate;
dopts.maxwidth = 2*mot.samplingRate;

dopts.frames2skip_begin = 0;
dopts.frames2skip_end   = 0;
dopts.allowedSteps      = [ 1 1; ...
                            1 0; ...
                            0 1; ...
                            1 2; ...
                            2 1; ...
                            2 2; ...
                            5 5; ...
                            10 10; ...
                            ];
                            
switch nargin
    case 2
        options = dopts;
    case 3
        options = mergeOptions(varargin{1},dopts);
    otherwise
        error('wrong num of args');
end

%COMMENT nearest neighbors subset size
options.k = min(mot.nframes, options.k);
options.k_s1     = options.k;
options.k_s3     = options.k;

options.radius = options.generalized_radius * sqrt( numel(options.frame_offsets) );
options.radius_activities = options.generalized_radius_activities * sqrt( numel(options.frame_offsets) );
options.radius_s1 = options.radius_activities;
options.radius_s3 = options.radius;

k = options.k;

%% extract feature set
tic;

fmat = tw_extractFeatureSetFromMot(skel, mot, options.feature_set, options.frame_offsets);

timings.Create_Feature_Set = toc;
fprintf('Create Feature Set\t\t\tcompleted in % 2f seconds.\n', timings.Create_Feature_Set);

%% knn search
tic;

handle = ann_buildTree(fmat);
%COMMENT nnidx: NN frame index for each frame
%COMMENT nndists: NN distance for each frame
%COMMENT the frames are compared to all other frames in the motion sequence
[nnidx,nndists] = ann_queryTree(handle,fmat,k,'search_sch','fr','radius',options.radius);
ann_cleanup();
nnidx = cast(nnidx,'double');
nnidx(nnidx==0)=nan;

timings.Find_Motion_Similarities = toc;
fprintf('Find Motion Similarities\tcompleted in % 2f seconds.\n', timings.Find_Motion_Similarities);
    

%% step 1
tic;

%COMMENT regard only a subset of NN
nnidx_s1 = nnidx(1:options.k_s1,:);
nndist_s1 = nndists(1:options.k_s1,:);
%COMMENT filter out NN with a distance larger than radius_s1
[ nnidx_s1, nndist_s1 ] = tw_filterRadius( options.radius_s1, nnidx_s1, nndist_s1 );

[nnidx_s1 nndist_s1 immat_s1] = ComputeWeightedNeighbourhood(fmat,nnidx_s1,nndist_s1,options.minrad,options.maxrad);

%COMMENT get forward and backward cuts
[~,cutsf,cutsb] = tw_countNNperSquare_old(nnidx_s1);

% Join fwd and bwd cuts:
% cuts = sortOutCuts([cutsf cutsb],2*options.minwidth);
%COMMENT join forward and backward cuts
cuts = sort([cutsf cutsb]);

timings.Find_Motion_Activities = toc;
fprintf('Find Motion Activities\t\tcompleted in % 2f seconds.\n', timings.Find_Motion_Activities);

%% Step 2
tic;

% Search fur subsegments:

% adjust neighbourhoods to respective velocity
[nnidx_s2 nndists_s2 immat_s2] = ComputeWeightedNeighbourhood(fmat,nnidx,nndists,options.minrad,options.maxrad);

% % % if options.DEBUG
% % %     figure();
% % %     colormap gray;
% % %     imagescnan([],[],immat2',isnan(immat2'),'b');
% % %     title(mot.filename,'interpreter','none');
% % % end

Afull = buildDAGmatrix(nnidx_s2,nndists_s2,options);

allsubcuts = [];
allpathlengths = [];
for i=1:numel(cuts)+1
    if i==1
       sf = 1;
    else
        sf = cuts(i-1);
    end
    if i==numel(cuts)+1
       ef = mot.nframes;
    else
       ef = cuts(i);
    end
  
%     curidx = curidx-sf+1;
    
    [Ai,curidx,curdis]  = cutDAGMatrixSquare(Afull,sf,ef,size(nnidx_s2,1),nnidx_s2,nndists_s2);
   
%     Ai = buildDAGmatrix(curidx,curdis,options);
    
    [paths,dists] = getSubsegments(curidx,curdis,Ai,options,sf);
    
    %check for soft symmetry: any valid path except for first diagonal
    %should have a pseudo-symmetrical twin (use first diagonal as symmtry axis)
    [paths, dists] = checkSoftSymmetry(paths,dists,sf);

    subcuts = zeros(1,numel(paths));
    subpaths = subcuts;
    for j=1:numel(paths)
        
        pcurlen  = paths{j}(1,1)-paths{j}(end,1)+1;
        pcomplen = paths{j}(1,2)-paths{j}(end,2);
        lval     = pcomplen/pcurlen;
        
        %check for very short paths
        if size(paths{j},1)>=options.minwidth && lval > 0.5 && lval < 2
            subcuts(j)=paths{j}(end,2);
            subpaths(j) = size(paths{j},1);%/pcomplen;
        end
        
    end
    [subcuts, ia, ic] = unique(subcuts);
    subpaths = subpaths(ia);
    allsubcuts = [allsubcuts subcuts];
    allpathlengths = [allpathlengths subpaths];
end

sortopts.priorities = cuts;
sortopts.paths = [zeros(1,numel(cuts)) allpathlengths];
allsubcuts = sortOutCuts([cuts allsubcuts],options.minwidth,sortopts);

allsubcuts = allsubcuts(2:end);
allpathlengths = allpathlengths(2:end);
sortopts.paths = [zeros(1,numel(cuts)) allpathlengths];

timings.Find_Motion_Segments = toc;
fprintf('Find Motion Segments\t\tcompleted in % 2f seconds.\n', timings.Find_Motion_Segments);


%% Step 3
tic;

% Cluster results of second step

c_options = options;
clmat  = zeros(numel(allsubcuts)+1);

%list of indices belonging to cuts in allsubcuts
foo = computeInDe(allsubcuts,mot.nframes);
% foo = foo(:,2:end); foo(2,1) = 1;
inbetweeners = computeInDe(cuts,mot.nframes);



%if cursub is on the list of 'inbetweeners' there should only be connecting edges
% linking cursub to other elements from the same list! 
if ~isequal(ef,mot.nframes)
    if ismember(sf,sortopts.priorities)&& ismember(ef+1,sortopts.priorities)
        allsubcutstmp = sortopts.priorities;
        footmp = inbetweeners;
    else
        allsubcutstmp = allsubcuts;
        footmp = foo;
    end
else
    allsubcutstmp = allsubcuts;
    footmp = foo;
end
allsubcutstmp = [allsubcutstmp mot.nframes];

nnidx_s3 = nnidx(1:options.k_s3,:);
nndists_s3 = nndists(1:options.k_s3,:);
nnidx_s3(nndists_s3>options.radius_s3)=nan;

[nnidx_s3 nndists_s3 immat_s3] = ComputeWeightedNeighbourhood(fmat,nnidx_s3,nndists_s3,options.minrad,options.maxrad);

Afull = buildDAGmatrix(nnidx_s3,nndists_s3,options);

for cursub = 1:numel(allsubcuts)
    %compute current index and relevant frames
    cur_idx = foo(2,cursub):foo(3,cursub);
    sf = foo(2,cursub);
    ef = foo(3,cursub);
   
    %consider strip of Afull corresponding to current segment and look for
    %suitable paths
    
    [A_strip, str_idx, str_dists] = cutDAGMatrixRect(Afull,sf,ef,1,mot.nframes,options.k,nnidx_s3,nndists_s3);
    [str_paths,~] = getSubsegments(str_idx,str_dists,A_strip,options);

% % %     if options.DEBUG
% % %         %draw self similarity matrix
% % %         immat = nan(size(str_idx,2),mot.nframes);
% % %         for i=1:size(str_idx,2)
% % %             immat(i,str_idx(~isnan(str_idx(:,i)),i))=str_dists(~isnan(str_idx(:,i)),i) ;
% % %         end
% % %         % % 
% % %         figure()
% % %         colormap gray;
% % %         imagescnan([],[],immat,isnan(immat),'b');
% % % 
% % %         plotPaths(str_paths,'red');
% % %         ylim = get(gca,'ylim');
% % %         for i=1:numel(allsubcuts)
% % %            line([allsubcuts(i) allsubcuts(i)], ylim,'color','white');
% % %         end
% % %     end 
    
    csids = getSegIds_local(cursub,allsubcutstmp);
    
    for c_path = 1:numel(str_paths)
        curpath=str_paths{c_path};
        len1 = curpath(1,1)-curpath(end,1);
        len2 = curpath(1,2)-curpath(end,1);
        c_clmat = computePathValidity(curpath,len1,len2);
        
        if c_clmat
            curids = curpath(end,1):curpath(1,1)+csids(1);
            refids = curpath(end,2):curpath(1,2);
            %find subcut where this path belongs
            [eef] = find(allsubcutstmp>=refids(end),1,'first'); % lest segment covered by path
            [ssf] = find(allsubcutstmp<=refids(1),1,'last')+1;    % first segment covered by path
            if isempty(ssf)
                ssf = 1;
            end
            for cref = ssf:eef
                rsids = getSegIds_local(cref,allsubcutstmp);
                if getcoverage_local(refids,curids,rsids,csids,0.5)
                    clmat(cursub,cref) = 1;
                end
            end
            
            
% %             if ~isequal(eef,1)
% %                 [ssf] = find(allsubcutstmp<=refids(1),1,'last');
% %                 clmat(cursub,ssf+1:eef-1) = 1;
% %             end

        end
    end

%     for compsub = 1:numel(allsubcutstmp)
%         
%         fprintf('\b\b\b\b\b\b%02i-%02i\n',cursub,compsub);        
%         sidx = footmp(2,compsub);
%         eidx = footmp(3,compsub);
%     
%         [Ai,curidx,curdis]  = cutDAGMatrixRect(Afull,sf,ef,sidx,eidx,options.k,nnidx_s2,nndists_s2);
%         [allpaths,alldists] = getSubsegments(curidx,curdis,Ai,options);
%         
%         allsizes = cell2mat(cellfun(@(x) size(x,1), allpaths, 'UniformOutput', false));
%         maxsize=find(allsizes==max(allsizes),1,'first');
%         
%         valpath = [allpaths{maxsize,1}];
%         cmp = foo(2,cursub):foo(3,cursub);
%         comp = footmp(2,compsub):footmp(3,compsub);
%         
%         curlen = numel(cmp);
%         complen = numel(comp);
%         
%         if ~isempty(valpath)
%             
%             pcurlen  = valpath(1,1)-valpath(end,1);
%             pcomplen = valpath(1,2)-valpath(end,2);
%             lval     = pcomplen/pcurlen;
%             xcover   = pcurlen /curlen;
%             ycover   = pcomplen/complen;
%             
%         else
%             lval = 0;
%             xcover = 0;
%             ycover = 0;
%         end
% 
%         pp = 0.30;
%         
%         sl = 1.3;
%         if ycover > pp && xcover > pp && ~(lval > sl || lval < 1/sl)
%             clmat(cursub,compsub) = 1;%alldists(maxsize);
%         end
%         
%         
% %         clmat(cursub,compsub) = sum(sum(~isnan(curcompidx)));
% %         clmat2(cursub,compsub) = clmat(cursub,compsub)/(abs(foo(2,cursub)-foo(3,cursub))*abs(foo(2,compsub)-foo(3,compsub)));
% 
%     end

end

clmat = sparse(clmat);
comps = components_mex(clmat);


%%
%DEBUG PLOTS

%create first plot representing the two segmentation steps 
if options.DEBUG
    %draw self similarity matrix
    immat = nan(size(fmat,2));
    for i=1:size(fmat,2)
        immat(i,nnidx_s1(~isnan(nnidx_s1(:,i)),i))=nndist_s1(~isnan(nnidx_s1(:,i)),i) ;
    end
    % % 
    figure()
    subplot(1,2,1)
    colormap gray;
    imagescnan([],[],immat_s1,isnan(immat_s1),'b');
    title(mot.filename,'interpreter','none');
    set(gcf,'position',[100 10 1800 900]);
    drawnow();
    
    %plot forward and backward plots separately
    
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    for i=1:numel(cuts)
        
        line([cuts(i) cuts(i)], ylim,'color','green','linewidth',2);
        line(xlim, [cuts(i) cuts(i)],'color','green','linewidth',2);
    end
    
    
    subplot(1,2,2)
    colormap gray;
    imagescnan([],[],immat_s2,isnan(immat_s2),'b');
    title(mot.filename,'interpreter','none');
    drawnow();
    
    for i=1:numel(allsubcuts)
        line([allsubcuts(i) allsubcuts(i)], ylim,'color','white');
        line(xlim, [allsubcuts(i) allsubcuts(i)],'color','white');
    end
    for i=1:numel(cuts)
        
        line([cuts(i) cuts(i)], ylim,'color','green','linewidth',2);
        line(xlim, [cuts(i) cuts(i)],'color','green','linewidth',2);
    end
    title(mot.filename,'interpreter','none');
    drawnow();
   
    
end

    
    
%% OTHER PLOTS


colors = jet(max(comps));

mots = cell(numel(cuts),1);
for i=1:numel(cuts)+1
   if i==1
       sf = 1;
   else
        sf = cuts(i-1);
   end
   
   if i==numel(cuts)+1
       ef = mot.nframes;
   else
       ef = cuts(i);
   end
   
   mots{i} = cutMotion(mot,sf,ef);
end


submots = cell(numel(allsubcuts),1);
for i=1:numel(allsubcuts)+1
   if i==1
       sf = 1;
   else
        sf = allsubcuts(i-1);
   end
   
   if i==numel(allsubcuts)+1
       ef = mot.nframes;
   else
       ef = allsubcuts(i);
   end
   
   submots{i} = cutMotion(mot,sf,ef);
   submots{i}.colors = colors(comps(i),:);
end

% % % if options.DEBUG
% % %     % create a new plot demonstrating clustering step (together with primitives)
% % %     figure;
% % % %     ax1=subplot(2,2,[1 3]);
% % %     colormap gray;
% % %     pos1 = get(ax1,'Position');
% % %     set(ax1,'Position',[pos1(1) pos1(2) pos1(3) pos1(4)]); 
% % %     imagescnan([],[],immat_s3,isnan(immat_s3),'b');
% % %     title(mot.filename,'interpreter','none');
% % %     ax3=subplot(2,2,4);
% % % %     pos2 = get(ax2,'Position');
% % % %     set(ax2,'Position',[pos2(1) pos2(2) pos2(3) 0.1*pos2(4)]);  
% % %     pos3 = get(ax3,'Position');
% % %     set(ax3,'Position',[pos3(1) pos3(2) pos3(3) 0.1*pos3(4)]);
% % %     plotCuts(ax3,[allsubcuts mot.nframes],[comps; comps(1)]);
% % %     ax2=subplot(2,2,2);
% % %     colormap gray;
% % %     pos2 = get(ax2,'Position');
% % %     set(ax2,'Position',[pos2(1) pos3(2)+0.2*pos3(4) pos2(3) 2.185*pos2(4)]); 
% % %     imagescnan([],[],immat_s3,isnan(immat_s3),'b');
% % %     title(mot.filename,'interpreter','none');
% % %     for i=1:numel(allsubcuts)
% % %         line([allsubcuts(i) allsubcuts(i)], ylim,'color','white');
% % %         %         line(xlim, [allsubcuts(i) allsubcuts(i)],'color','white');
% % %     end
% % %     for i=1:numel(cuts)
% % %         
% % %         line([cuts(i) cuts(i)], ylim,'color','green','linewidth',2);
% % %         line(xlim, [cuts(i) cuts(i)],'color','green','linewidth',2);
% % %     end
% % %     drawnow()
% % % end


timings.Cluster_Motion_Segments = toc;
fprintf('Cluster Motion Segments\t\tcompleted in % 2f seconds.\n', timings.Cluster_Motion_Segments);

end

%% HELPER FUNCTIONS (tidy up!!)

function g = pathgradient(p)

    g = abs(p(1,1)-p(end,1)+1)/abs(p(1,2)-p(end,2)+1);

end

function segs = convertPathsToSegments_local(allpaths,alldists)
if ~isempty(allpaths)
    nframes = 1;
    for i=1:numel(allpaths)
        if nframes < allpaths{i}(end,1)
            nframes = allpaths{i}(end,1);
        end
    end
    
    signal = zeros(1,nframes+2);
    for i=1:numel(allpaths)
        signal(allpaths{i}(end,1)+1) = alldists(i);
    end
    
    [~,segs] = findpeaks(signal,'MINPEAKDISTANCE',min(15,numel(signal)-1));
    segs = segs-1;
else
    segs = [];
end
end

function plotrowhists_local(immat)

nbuckets = 16;
immat(isnan(immat))=-1;

hists = zeros(nbuckets,size(immat,2));

for i=1:size(immat,2)
   hists(:,i) = hist(immat(:,i),nbuckets); 
end

figure();
imagesc(log(hists));

end



function plotMultiResMats(immat)
figure()
nLevels = 10;

nanrep = max(immat(:))*2;
immat(isnan(immat))=nanrep;

levels = 2.^(0:nLevels-1);
lcount = 1;
for curLevel=levels
   
   width = floor(size(immat,2)/curLevel);
    
   curM = zeros(curLevel);

   for r=1:curLevel
       for c=1:curLevel
           
           tmp = (immat( (r-1)*width+1:r*width,(c-1)*width+1:c*width ));
           curM(r,c) = var(tmp(:));
           %curM(r,c) = mean(mean(immat( (r-1)*width+1:r*width,(c-1)*width+1:c*width )));

       end
   end
  
   rows = ceil(nLevels/5);
   cols = min(5,nLevels);
   
   subplot(rows,cols,lcount);
   imagesc(curM);
   set(gca,'clim',[0 nanrep]);
   lcount = lcount+1;

    
end

colorbar();

end

function [immat] = paintExcludedBlocks(mot,immat,blockgraph)



for lev = 1:blockgraph.numLevels
   clevel = blockgraph.getLvlNodes(lev);
   for cn = clevel
       r = blockgraph.nodes{cn}.nodeMatIDsR;
       c = blockgraph.nodes{cn}.nodeMatIDsC;
       if isequal(blockgraph.nodes{cn}.dealtwith,1)
           immat(r,c) = 0;
       else
           immat(r,c) = 1;
       end
           
   end
end

fh=figure();
colormap gray
imagescnan([],[],immat,isnan(immat),[0.2 0 0.8]);
title([mot.filename ' at ' num2str(mot.samplingRate) 'Hz'],'interpreter','none');
xlabel('time [frames]')
ylabel('time [frames]')

drawnow();

end

function res = sortOutCuts(c,p,opts)

[res ia] = unique(c);

pathlengths = zeros(1,numel(c));
priorities = 0;

if nargin == 3
    if isfield(opts,'priorities')
        priorities = opts.priorities;
    end
    
    if isfield(opts,'paths')
        pathlengths = opts.paths;
        pathlengths = pathlengths(ia);
    end
end

i=1;
n = size(res,2);
while i < n-1
    dr=abs(res(i)-res(i+1));
    curme = [res(i) res(i+1)];
    if res(i)==0
        res(i) =  [];
        pathlengths(i) = [];
        n=n-1;
    elseif dr < p
        if numel(intersect(curme,priorities)) == 2
            n=n+1;
        elseif ismember(res(i),priorities)
            res(i+1)=[];
            pathlengths(i+1)=[];
        elseif ismember(res(i+1),priorities)
            res(i)=[];
            pathlengths(i) =[];
        elseif pathlengths(i)> pathlengths(i+1)
            res(i+1)=[];
            pathlengths(i+1) =[];
        else
            res(i) = [];
            pathlengths(i) =[];
        end
        n=n-1;
    end
    
    if abs(res(i)-res(i+1)) >= p || numel(intersect(curme,priorities)) == 2
        i = i+1;
    end

end
end

function [InDe]=computeInDe(cuts,n)

InDe = zeros(3,numel(cuts));

for cursub = 1:numel(cuts)
    if isequal(cursub, numel(cuts))
        ef = n;
    else
        ef = cuts(cursub)-1;
    end
    if isequal(cursub, 1)
        sf = 1;
    else
        sf = cuts(cursub-1);
    end
    InDe(:,cursub) = [cuts(cursub); sf; ef];
end

end


function coverage = getcoverage_local(pr,pc,sr,sc,P)

coverage = 1; % Asssume path is ok in the first place.

xcov = sum(ismember(sr,pr))/numel(sr);
ycov = sum(ismember(sc,pc))/numel(sc);

if xcov<P
    coverage = 0;
elseif ycov<P
    coverage = 0;
end

end

function ids = getSegIds_local(c,cuts)

if c==1
    sf = 1;
else
    sf = cuts(c-1);
end

ef = cuts(c);

ids = sf:ef;

end

