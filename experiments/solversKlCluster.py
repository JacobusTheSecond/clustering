from solvers import DriftersSolver
import klcluster as kl
import numpy as np

class KlClusterDriftersSolver(DriftersSolver):
    def __init__(self, drifterFiles):
        super().__init__(drifterFiles)

        self.curves = kl.Curves()
        for i, curvedata in enumerate(self.datacurves):
            self.curves.add(kl.Curve(curvedata, f"Curve_{i}"))

        self.DELTA = 15000
        self.COMPLEXITY = 10
        self.ROUNDS = 1

        self.cc = kl.CurveClusterer()
        self.cc.initCurves(self.curves, self.DELTA)
        self.simplifiedCurves = self.cc.getSimplifications()

    def solve(self):
        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS)

        # convert cluster centers to np array
        self.clustercurves = []
        for cluster in self.clusters:
            curve = []

            center = cluster.center()
            start = int(center.start.value)
            end = int(center.end.value)
            curve = self.simplifiedCurves[center.curve]

            curvedata = []
            for i in range(start, end+1):
                if (i < 0 or i >= len(curve)):
                    raise Exception("Index out of bounds")
                # print(curve[i])
                curvedata.append(curve[i].values)

            self.clustercurves.append(np.array(curvedata))
