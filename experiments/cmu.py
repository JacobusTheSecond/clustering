from cmuSolvers import KlClusterCMUSolver

solver = KlClusterCMUSolver(1)
segments, segmentation = solver.solve()

solver.plotSegmentsAndGT(segments)