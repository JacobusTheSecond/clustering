% iterativeForwKinematics is faster than the recursive version!!

function [jointTrajectories,localSystems] = iterativeForwKinematics(skel,mot)

%%% Zeilen einkommentieren um Fehler abzufangen,
%%% auskommentieren um Geschwindigkeit zu steigern
% if isempty(mot.rootTranslation)
%     mot.rootTranslation = zeros(3,mot.nframes);
% end
% 
% if isempty(mot.rotationQuat) || isempty(mot.rotationQuat{1})
%     if ~isempty(mot.rotationEuler) || ~isempty(mot.rotationEuler{1})
%         mot = convert2quat(skel,mot);
%     else
%         error('Rotation data unavailable.');
%     end
% end
    
if ~iscell(mot.rotationQuat)
    tmp                    = mat2cell( mot.rotationQuat, ...
        4*ones(1,size(mot.rotationQuat,1)/4), ...
        size(mot.rotationQuat,2));
    mot.rotationQuat = cell(skel.njoints,1);
    mot.rotationQuat(skel.animated) = tmp;
end



localSystems         = cell(skel.njoints,1);
localSystems{1}      = mot.rotationQuat{1};
jointTrajectories    = cell(skel.njoints,1);
jointTrajectories{1} = mot.rootTranslation;

for i=1:size(skel.paths,1)
    for j=2:numel(skel.paths{i})
        joint = skel.paths{i}(j);
        pred  = skel.paths{i}(j-1);
        
        if isempty(mot.rotationQuat{joint})
            localSystems{joint,1} = localSystems{pred};
        else
            localSystems{joint,1} = C_quatmult(localSystems{pred},mot.rotationQuat{joint});
        end
        jointTrajectories{joint} = jointTrajectories{pred}...
                                     + C_quatrot(skel.nodes(joint).offset,localSystems{joint});
    end
end