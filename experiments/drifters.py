import os

from drifterSolvers import KlClusterDriftersSolver

klcluster = True
datafolder = os.path.join(os.path.dirname(__file__), "../data_drifters/world3d_txt")
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]

print(drifterFiles[9:10])

solver = KlClusterDriftersSolver(drifterFiles[9:10]) if klcluster else None

result = solver.solve()
solver.plotResults()