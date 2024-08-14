import os

from drifterSolvers import KlClusterDriftersSolver

n = 1000

klcluster = True
datafolder = os.path.join(os.path.dirname(__file__), "../data_drifters/world3d_txt")
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]

solver = KlClusterDriftersSolver(drifterFiles[0:n]) if klcluster else None

result = solver.solve(onlyRelevantClusters=True)
solver.plotResults()
