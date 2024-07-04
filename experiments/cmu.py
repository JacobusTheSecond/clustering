from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver

TAG = 2
klClusterSolver = KlClusterCMUSolver(TAG)
tmmCMUSolver = TmmCMUSolver(TAG)
gmmSolver = AcaCMUSolver(TAG, "gmm")
acaSolver = AcaCMUSolver(TAG, "aca")
hacaSolver = AcaCMUSolver(TAG, "haca")


segmentsKlCluster = klClusterSolver.solve()
segmentsTmm = tmmCMUSolver.solve()
#print(klClusterSolver.calculateAccuracy(segmentsKlCluster, TAG))

#klClusterSolver.plotSegmentsAndMatching()

segmentsGmm = gmmSolver.solve()

segmentsAca = acaSolver.solve()
segmentsHaca = hacaSolver.solve()

klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsTmm, segmentsAca, segmentsHaca, segmentsGmm], ["KlCluster", "Tmm", "Aca", "Haca", "Gmm"])
