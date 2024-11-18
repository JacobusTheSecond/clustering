from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver

TAG = 1
klClusterSolver5 = KlClusterCMUSolver(TAG, 0.8, 1.35, 5)
klClusterSolver10 = KlClusterCMUSolver(TAG, 0.8, 1.7, 10)
klClusterSolver15 = KlClusterCMUSolver(TAG, 0.8, 1.05, 15)
tmmCMUSolver = TmmCMUSolver(TAG)
gmmSolver = AcaCMUSolver(TAG, "gmm")
acaSolver = AcaCMUSolver(TAG, "aca")
hacaSolver = AcaCMUSolver(TAG, "haca")


segmentsKlCluster5 = klClusterSolver5.solve(mergeOverlappingClusters=False)
segmentsKlCluster5 = klClusterSolver5.solveAndTime(timing_rounds=3, mergeOverlappingClusters=False)
#segmentsKlCluster10 = klClusterSolver10.solve(mergeOverlappingClusters=False)
#segmentsKlCluster15 = klClusterSolver15.solve(mergeOverlappingClusters=False)
segmentsTmm = tmmCMUSolver.solve()

acc = tmmCMUSolver.calculateAccuracy(segmentsTmm, TAG)

print(acc, "acc")
#segmentsGmm = gmmSolver.solve()
#segmentsAca = acaSolver.solve()
#segmentsHaca = hacaSolver.solve()

#klClusterSolver5.plotSegmentsAndGT([segmentsKlCluster5, segmentsKlCluster10, segmentsKlCluster15, segmentsTmm, segmentsAca, segmentsHaca, segmentsGmm], ["KlCluster5", "KlCluster10", "KlCluster15", "Tmm", "Aca", "Haca", "Gmm"])
