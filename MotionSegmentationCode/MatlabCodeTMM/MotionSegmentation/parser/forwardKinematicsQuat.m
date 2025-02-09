function jointTrajectories = forwardKinematicsQuat(skel,mot)
jointTrajectories = recursive_forwardKinematicsQuat( ...
                     skel,...
                     mot,...
                     1,...
                     mot.rootTranslation + repmat(skel.nodes(1).offset,1,mot.nframes),...
                     quatmult( ...
                               repmat(skel.rootRotationalOffsetQuat,1,mot.nframes), ...
                               mot.rotationQuat{1}),...
                     mot.jointTrajectories);