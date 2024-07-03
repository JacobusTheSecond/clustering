% iterativeForwKinematics is faster than the recursive version!!

function [jointTrajectories,localSystems] = iterativeForwKinematics_mats(skel,rootTranslation,rotationQuat)

nrOfFrames = size(rootTranslation,2);

localSystems             = zeros(skel.njoints*4,nrOfFrames);
localSystems(1:4,:)      = rotationQuat(1:4,:);
jointTrajectories        = zeros(skel.njoints*3,nrOfFrames);
jointTrajectories(1:3,:) = rootTranslation;

if size(rotationQuat,1)==numel(skel.animated)*4
    id = [1;0;0;0];
    for i=sort(skel.unanimated)'
        l = rotationQuat(1:(i-1)*4,:);
        u = rotationQuat((i-1)*4+1:end,:);
        rotationQuat = [l;id(:,ones(1,nrOfFrames));u];
    end
end

for i=1:size(skel.paths,1)
    for j=2:numel(skel.paths{i})
        joint = skel.paths{i}(j);
        pred  = skel.paths{i}(j-1);
        localSystems(joint*4-3:joint*4,:) = C_quatmult(localSystems(pred*4-3:pred*4,:),rotationQuat(joint*4-3:joint*4,:));
        jointTrajectories(joint*3-2:joint*3,:) = jointTrajectories(pred*3-2:pred*3,:)...
                                     + C_quatrot(skel.nodes(joint).offset,localSystems(joint*4-3:joint*4,:));
    end
end