from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver

klClusterSolver = KlClusterCMUSolver(1)
gmmSolver = AcaCMUSolver(1, "gmm")
acaSolver = AcaCMUSolver(1, "aca")
hacaSolver = AcaCMUSolver(1, "haca")

segmentsKlCluster = klClusterSolver.solve()
segmentsGmm = gmmSolver.solve()
segmentsAca = acaSolver.solve()
segmentsHaca = hacaSolver.solve()

klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsGmm, segmentsAca, segmentsHaca], ["KlCluster", "Gmm", "Aca", "Haca"])