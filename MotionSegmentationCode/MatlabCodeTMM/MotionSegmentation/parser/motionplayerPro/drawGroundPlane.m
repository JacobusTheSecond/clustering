    function drawGroundPlane(varargin)
    
    global SCENE;
    
        if SCENE.status.groundPlane_drawn
            SCENE.status.groundPlane_drawn = false;
%             set(SCENE.handles.groundPlane,'Visible','off');
            set(SCENE.handles.drawGroundPlane_Button,'CData',SCENE.buttons.groundPlane+0.5,'TooltipString','draw ground plane');
            set(SCENE.handles.groundPlane,'FaceAlpha',0,'EdgeColor','none');
        else
            SCENE.status.groundPlane_drawn = true;
%             set(SCENE.handles.groundPlane,'Visible','on');
            set(SCENE.handles.drawGroundPlane_Button,'CData',SCENE.buttons.groundPlane,'TooltipString','hide ground plane');
            set(SCENE.handles.groundPlane,'FaceAlpha',0.7,'EdgeColor','none');
        end
    end