import os

from drifterSolvers import KlClusterDriftersSolver

n = 150

klcluster = True
datafolder = os.path.join(os.path.dirname(__file__), "../data_drifters/world3d_txt")
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]

solver = KlClusterDriftersSolver(drifterFiles[0:n], rossbyRadius=False) if klcluster else None

result = solver.solve(onlyRelevantClusters=True, withShow=False)
solver.drawInputField()
solver.drawResultField()
solver.plotInputAndResult()
solver.plotCurve()
solver.plotResults()
