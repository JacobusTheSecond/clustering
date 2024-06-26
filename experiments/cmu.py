from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver

TAG = 1
klClusterSolver = KlClusterCMUSolver(TAG)
gmmSolver = AcaCMUSolver(TAG, "gmm")
acaSolver = AcaCMUSolver(TAG, "aca")
hacaSolver = AcaCMUSolver(TAG, "haca")


segmentsKlCluster = klClusterSolver.solve()
#print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG))

segmentsGmm = gmmSolver.solve()

segmentsAca = acaSolver.solve()
segmentsHaca = hacaSolver.solve()

klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsGmm, segmentsAca, segmentsHaca], ["KlCluster", "Gmm", "Aca", "Haca"])
