from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver

TAG = 4


klClusterSolver = KlClusterCMUSolver(TAG)
tmmCMUSolver = TmmCMUSolver(TAG)
#gmmSolver = AcaCMUSolver(TAG, "gmm")
#acaSolver = AcaCMUSolver(TAG, "aca")
#hacaSolver = AcaCMUSolver(i, "haca")
segmentsKlCluster = klClusterSolver.solve()
segmentsTmm = tmmCMUSolver.solve()


#segmentsTmm = tmmCMUSolver.solve()
#print(klClusterSolver.calculateAccuracy(segmentsKlCluster, TAG))

#klClusterSolver.plotSegmentsAndMatching()

#segmentsGmm = gmmSolver.solve()

#segmentsAca = acaSolver.solve()
#segmentsHaca = hacaSolver.solve()

klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsTmm], ["KlCluster", "tmm"])#, segmentsTmm, segmentsAca, segmentsHaca, segmentsGmm], ["KlCluster", "Tmm", "Aca", "Haca", "Gmm"])
#klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsTmm, segmentsAca, segmentsHaca, segmentsGmm], ["KlCluster", "Tmm", "Aca", "Haca", "Gmm"])
