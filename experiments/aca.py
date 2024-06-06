import matlab.engine
import os

eng = matlab.engine.start_matlab()  #connect_matlab()
eng.cd(f'{os.getcwd()}/../aca/aca', nargout=0)
# eng.make(nargout=0)
eng.addPath(nargout=0) # add aca paths like src or lib
eng.addpath("..") # add runAca, runGmm, ... to path
gt, seg = eng.runGmm(1, nargout=2)
print(gt)
print(seg)
## stop matlab.engine
eng.quit()

