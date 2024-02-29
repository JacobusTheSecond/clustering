import numpy as np
import klcluster as kl
cs = kl.Curves()
for i in range(10):
	cs.add(kl.Curve(np.random.rand(10,3)*100))
cc = kl.CurveClusterer()
cc.initCurves(cs,1)
#cr1 = cc.greedyCover(5,1)
cr2 = cc.greedyCoverAgressiveFilter(5,1)
#print(cr1.__len__(),cr2.__len__())
for c in cr2:
	print("Center:",c.center().values())
	print("Matching:")
	print(c.values())
	print("----")
