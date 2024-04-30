import numpy as np
import klcluster as kl

trials = [np.genfromtxt("data/86_"+str(i)+".txt",delimiter=" ") for i in range(1,2)]
cs=kl.Curves()
for trial in trials:
	cs.add(kl.Curve(trial))

cc=kl.CurveClusterer()
cc.initCurves(cs,1.25)
cr = cc.greedyCover(10,1)
for cluster in cr:
	print("Center:",cluster.center().values())
	print("Matching:")
	print(cluster.values())
	print("------")
