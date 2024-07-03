function playFunction(varargin)

global SCENE;

if (SCENE.status.curFrame == SCENE.nframes)
    SCENE.status.curFrame = 1;
end
SCENE.status.reverse    = false;
SCENE.status.running    = true;
SCENE.status.timeStamp  = (SCENE.status.curFrame)/SCENE.samplingRate;
SCENE.timeOffset        = SCENE.status.timeStamp;

tic;
while SCENE.status.running && (SCENE.status.curFrame<SCENE.nframes || SCENE.status.looped)
    nextFrame               = getNextFrame_local();
    if SCENE.status.running
        setFramePro(nextFrame);
        drawnow;
        SCENE.status.timeStamp = SCENE.timeOffset+toc*SCENE.status.speed;
    end
end
SCENE.status.running = false;
end

function nextFrame = getNextFrame_local()
    global SCENE;
    
    framesToDrop = SCENE.status.timeStamp*SCENE.samplingRate-SCENE.status.curFrame;
    
    if SCENE.status.reverse
        nextFrame = SCENE.status.curFrame+round(framesToDrop)-1;
    else
        nextFrame = SCENE.status.curFrame+round(framesToDrop)+1;
    end

    if nextFrame<=0
        if SCENE.status.looped
            nextFrame   = SCENE.nframes;
        else
            nextFrame = 1;
        end
    elseif nextFrame>SCENE.nframes
        if SCENE.status.looped
            nextFrame   = 1;%mod(nextFrame,SCENE.nframes);
        else
            nextFrame = SCENE.nframes;
        end
    end

    if (framesToDrop<0 && ~SCENE.status.reverse) || (framesToDrop>0 && SCENE.status.reverse)
        pause(abs(framesToDrop/SCENE.samplingRate));
    end
end