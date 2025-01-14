import os

from drifterSolvers import KlClusterDriftersSolver

#n = 1000
ns = [25,50,100,200,400,800,1600,3200]

klcluster = True
datafolder = os.path.join(os.path.dirname(__file__), "../data_drifters/world3d_txt")
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]
lowerbounds = []
solutionsizes = []

for n in ns:
    solver = KlClusterDriftersSolver(drifterFiles[0:n], rossbyRadius=False) if klcluster else None
    lowerbounds.append(solver.lowerbound(False))
    solutionsizes.append(solver.solutionSize(False))
    print(f"ns: {ns}")
    print(f"Lowerbounds: {lowerbounds}")
    print(f"Size of Sol: {solutionsizes}")
    #result = solver.solve(onlyRelevantClusters=True, withShow=False)
#solver.drawInputField()
#solver.drawResultField()
#solver.plotInputAndResult()
#solver.plotCurve()
#solver.plotResults()
