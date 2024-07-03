function SCENE = defaultSCENE()

SCENE.position      = [150,150];
% SCENE.size          = [480, 640]/2;
% SCENE.size          = [640 ,480];
% SCENE.size          = [1920/2 ,1080]; % viertel HD
% SCENE.size          = [1920 ,1080]; % Full HD
SCENE.size          = [1280 , 720];% 720p 1280�720
% SCENE.size          = [426, 720];
% SCENE.size          = [1920/3 ,1000]; % Spezial f�r 4 nebeneinander in figure*

SCENE.opened        =  fix(clock);

SCENE.nmots         = 0;
SCENE.nskels        = 0;
SCENE.nframes       = 0;
SCENE.boundingBox   = [-1;1;-1;1;-1;1]*150; % minimum bounding box
SCENE.samplingRate  = 0;
SCENE.objects       = [];

SCENE.status = struct(...
    'running',                  false,...
    'curFrame',                 1,...
    'reverse',                  false,...
    'speed',                    1,...
    'looped',                   false,...
    'spread',                   false,...
    'timeStamp',                0,...
    'timeOffset',               0,...
    'mainAxis',                 'y',...
    'groundPlane_drawn',        false,...
    'coordSyst_drawn',          false,...
    'localCoordSyst1_drawn',    false,...
    'localCoordSyst2_drawn',    false,...
    'sensorCoordSyst_drawn',    false,...
    'jointIDs_drawn',           false,...
    'light',                    'on'); % 'on' or 'off'

SCENE.colors = struct(...
    'singleSkel_FaceVertexData',            [0 0 1;0 1 0;0 0 1;0 1 0;0 0 1;1 0 0],...
    'singleSkel_FaceColor',                 [0 0 1],...[1 0.4 0.2],...'interp',...
    'singleSkel_FaceAlpha',                 0.7,...1,...0.2,...
    'singleSkel_EdgeColor',                 'none',...'interp',...
    'multipleSkels_FaceVertexData_start',   rgb2hsv([1 0 0]),... rgb2hsv([54/255 54/255 193/255]),...
    'multipleSkels_FaceVertexData_end',     rgb2hsv([0 0 1]),... 
    'multipleSkels_FaceColor',              'flat',...
    'multipleSkels_FaceAlpha',              0.8,...
    'multipleSkels_EdgeColor',              'none',...
    'points_start',                         rgb2hsv([0.0 0.8 0.3]),...rgb2hsv([0.0 1.0 0.0]),... % green
    'points_end',                           rgb2hsv([1.0 0.0 0.0]),... % red
    'buttons_group1',                       [0.5,0.8,0.95],...
    'buttons_group2',                       [0.5,0.8,0.95],...%[0.8,.2,0.8],...
    'groundPlane',                          [1 1 1;0.95 0.95 0.95],...
    'groundPlane_FaceAlpha',                0.7,...
    'backgroundColor',                      [1 1 1]);

SCENE.bones                         = 'sticks'; % allowed: diamonds, sticks, tubes, spheres
SCENE.spreadOffset                  = 150;
SCENE.ahead                         = [0;0;1];
SCENE.lengthOfLocalCoordSystAxes    = 10;
SCENE.scaleFactorForC3D             = 1;
SCENE.spreadFormation               = 'line'; % allowed: 'line', 'square', not allowed: 'turtle'
SCENE.lightPosition                 = [0 300 0];
SCENE.gravity                       = 9.812;
SCENE.groundPlaneSquareSize         = 50; % in cm!!!

SCENE.virtualSensors = defaultSensors();

SCENE.obj_default = struct(...
    'size',         1,...
    'color',        'red');


SCENE.handles = struct(...
    'fig',-1,...
    'groundPlane',-1,...
    'control_Panel',-1,...
    'goto_First_Button',-1,...
    'play_reverse_Button',-1,...
    'pause_Button',-1,...
    'play_Button',-1,...
    'goto_Last_Button',-1,...
    'slower_Button',-1,...
    'faster_Button',-1,...
    'quit_Button',-1,...
    'loop_Button',-1,...
    'axis_x_Button',-1,...
    'axis_y_Button',-1,...
    'axis_z_Button',-1,...
    'spread_Button',-1,...
    'render_Button',-1,...
    'status_Panel',-1,...
    'curFrameLabel',-1,...
    'curSpeedLabel',-1,...
    'help_Button',-1,...
    'sliderHandle',-1,...
    'drawLocalCoordSyst_Button',-1,...
    'drawSensorCoordSyst_Button',-1,...
    'coordX',-1,...
    'coordY',-1,...
    'coordZ',-1,...
    'light',-1);

SCENE.buttons = MPPbuttons();

SCENE.rendering = struct(...
    'filename',         'timestamp',... % timestamp
    'fps',              25,...%'full',... % full reads sampling rate from SCENE
    'compression',      'none',...
    'quality',          100,...
    'bgColor',          [235, 242, 242]);

end