from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver

TAG = 1
#klClusterSolver = KlClusterCMUSolver(TAG)
#gmmSolver = AcaCMUSolver(TAG, "gmm")
#acaSolver = AcaCMUSolver(TAG, "aca")
#hacaSolver = AcaCMUSolver(TAG, "haca")
tmmSolver = TmmCMUSolver(TAG)


#segmentsKlCluster = klClusterSolver.solve()
#print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG))

#segmentsGmm = gmmSolver.solve()

#segmentsAca = acaSolver.solve()
#segmentsHaca = hacaSolver.solve()
segmentsTmm = tmmSolver.solve()

print(segmentsTmm)

tmmSolver.plotSegmentsAndGT([segmentsTmm], ["Tmm"])
#klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsGmm, segmentsAca, segmentsHaca, segmentsTmm], ["KlCluster", "Gmm", "Aca", "Haca", "Tmm"])
