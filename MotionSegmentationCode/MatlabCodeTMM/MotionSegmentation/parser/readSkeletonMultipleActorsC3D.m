function [skel,mots] = readSkeletonMultipleActorsC3D(filename)

[skel,mot3d] = readMocapSmartLCS(filename,[],true);

mots = buildMotsFromNameMap_local(mot3d);

end

function mots3d = buildMotsFromNameMap_local(mot3d)

% Hard coded marker set HDM DB!
markers = { 'LFHD';'RFHD';'LBHD';'RBHD'; ...
            'C7';  'T10'; 'CLAV';'STRN'; ...
            'RBAC';'LBAC';'LSHO';'LUPA'; ...
            'LELB';'LFRM';'LWRA';'LWRB'; ...
            'LFIN';'RSHO';'RUPA';'RELB'; ...
            'RFRM';'RWRA';'RWRB';'RFIN'; ...
            'LFWT';'RFWT';'LMWT';'RMWT'; ...
            'LBWT';'RBWT';'LTHI';'LKNE'; ...
            'LSHN';'LANK';'LHEE';'LMT1'; ...
            'LMT5';'LTOE';'RTHI';'RKNE'; ...
            'RSHN';'RANK';'RHEE';'RMT1'; ...
            'RMT5';'RTOE'};

% Find all markers with these names:
idx = cell(numel(markers),1);

for curM = 1:numel(markers)
   
    idx{curM} = find(strncmp(markers(curM), mot3d.nameMap(:,1),numel(markers{curM})));
    
end

% check if all markers are there same often.
nmarkers = size(idx{1},1);
for i=2:numel(markers)
   if size(idx{i},1)~=nmarkers
       error('Not the same number of standard markers were found.')
   end
end

% now build mots
mots3d    = cell(nmarkers,1);
markerids = zeros(numel(markers),1);
for cm = 1:nmarkers
   
    % copy marker ids
    for i=1:numel(markers)
        markerids(i) = idx{i}(cm);
    end

    % copy mot
    mots3d{cm} = mot3d;
    
    % remove data from other actors
    mots3d{cm}.jointTrajectories = mots3d{cm}.jointTrajectories(markerids);
    mots3d{cm}.nameMap =  mots3d{cm}.nameMap(markerids,:);
    mots3d{cm}.njoints = numel(markers);
    mots3d{cm}.filename = [mots3d{cm}.filename '_' num2str(cm)];
end

end