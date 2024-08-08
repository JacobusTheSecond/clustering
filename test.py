import numpy as np
import klcluster as kl
cs = kl.Curves()

labelNum = 5

gts = kl.GroundTruths()

curveLen = 30

for i in range(10):
	c = kl.Curve(np.random.rand(curveLen,3)*10)
	w = c.getWeights()
	w*=3
	c.setWeights(w)
	cs.add(c)
	gtBase = np.sort(np.random.choice(range(curveLen),5,False))
	gtBase += curveLen-gtBase[-1]	
		

	gt = kl.GroundTruth()
	for j in gtBase:
		label = np.random.randint(0,labelNum)
		gt.add(label,j)
		#print(label,j)
	gts.add(gt)
#w = cs[0].getWeights()
#w/=2
#cs[0].setWeights(w)
#print(cs[0].getWeights())
#print(w)
#print(cs[0].getWeights() == w)
#print(len(cs))
cc = kl.CurveClusterer()
cc.initCurvesWithGT(cs,1,gts)

sgts = cc.getSimplifiedGTs()
for j in range(10):
	for i in range(5):
		print(sgts[j].labelAt(i),sgts[j].paramAt(i).value)

print(cc.getSimplifications())
cr1 = cc.greedyCover(5,1)
#cr2 = cc.greedyCoverAgressiveFilter(5,1)
#print(cr1.__len__(),cr2.__len__())
for c in cr1:
	print("Center:",c.center().values())
	print("Matching:")
	print(c.values())
	print("----")
