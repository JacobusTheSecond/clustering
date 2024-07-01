function res=convertMot2DBmat(skel,mot)

    res.nrOfFrames   = mot.nframes;
    res.frameRate    = mot.samplingRate;
    res.motStartIDs  = 1;
    res.origRootPos  = mot.rootTranslation;
    
%     if skel.njoints<31
%         mot_tmp = padSimpleMot(skel,mot);
%     else
        mot_tmp = mot;
%     end

    if ~isempty(mot_tmp.jointTrajectories)
        if iscell(mot_tmp.jointTrajectories)
            res.posOrig      = cell2mat(mot_tmp.jointTrajectories);
        else
            res.posOrig      = mot_tmp.jointTrajectories;
        end
        
        if mot.nframes>5
            mot_tmp = addVelToMot(mot_tmp);
            mot_tmp = addAccToMot(mot_tmp);
            res.velOrig = cell2mat(mot_tmp.jointVelocities);
            res.accOrig = cell2mat(mot_tmp.jointAccelerations);
        end
    end
    
    if isfield(mot,'jointVelocities')
        if iscell(mot.jointVelocities)
            res.velOrig = cell2mat(mot.jointVelocities);
        else
            res.velOrig = mot.jointVelocities;
        end
    end
    
    if isfield(mot,'jointAccelerations')
        if iscell(mot.jointAccelerations)
            res.accOrig = cell2mat(mot.jointAccelerations);
        else
            res.accOrig = mot.jointAccelerations;
        end
    end
    
    if iscell(mot.rotationQuat)
        rootRot = mot.rotationQuat{1};
    else
        rootRot = mot.rotationQuat(1:4,:);
    end

    if ~isempty(mot.rotationQuat)
        if size(rootRot,2)>1
            res.deltaRootPos = [ [0;0;0] C_quatrot(mot.rootTranslation(:,2:end)-mot.rootTranslation(:,1:end-1),C_quatinv(rootRot(:,2:end)))];
            res.deltaRootOri = [ [1;0;0;0] C_quatmult(rootRot(:,2:end),C_quatinv(rootRot(:,1:end-1)))];
            res.origRootOri  = rootRot;
        else
            res.deltaRootPos = [0;0;0];
            res.deltaRootOri = [1;0;0;0];
            res.origRootOri  = rootRot;
        end
    end

    mot.rootTranslation(:,:) = 0;
    
    if ~isempty(mot.rotationQuat)
        
        [mot,qy] = C_fitRootOrientationsFrameWise(skel,mot);
        res.invRootRot = qy;
        if isfield(mot,'rotationEuler')
            if isempty(mot.rotationEuler)
                mot = C_convert2euler(skel,mot);
            end
        end
        
%         if skel.njoints<31
%             mot = padSimpleMot(skel,mot);
%         end

        if iscell(mot.rotationQuat)
            res.quat = cell2mat(mot.rotationQuat(mot.animated));
        else
            res.quat = mot.rotationQuat;
        end
        
        if isfield(mot,'rotationEuler')
            res.euler = cell2mat(mot.rotationEuler(mot.animated));
        end
        
        res.pos = cell2mat(mot.jointTrajectories);
        if mot.nframes>5
            mot = addVelToMot(mot);
            mot = addAccToMot(mot);
            res.vel = cell2mat(mot.jointVelocities);
            res.acc = cell2mat(mot.jointAccelerations);
            res.velPlus = cell2mat(cellfun(@(x) C_quatrot(double(x),qy),mot.jointVelocities,'UniformOutput',0));
        else
            res.vel=zeros(mot.njoints*3,mot.nframes);
            res.acc=zeros(mot.njoints*3,mot.nframes);
        end
    end
    
    if isfield(mot,'phy');
        phycell = struct2cell(mot.phy);
        phymat  = cell2mat(phycell(2:end));
        res.phy = phymat;
    end
end