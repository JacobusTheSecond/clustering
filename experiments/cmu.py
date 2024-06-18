from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
TAG = 2
klClusterSolver = KlClusterCMUSolver(TAG)
gmmSolver = AcaCMUSolver(TAG, "gmm")
acaSolver = AcaCMUSolver(TAG, "aca")
hacaSolver = AcaCMUSolver(TAG, "haca")


segmentsKlCluster = klClusterSolver.solve()
print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG))

print("P")
segmentsGmm = gmmSolver.solve()
print("P")

segmentsAca = acaSolver.solve()
segmentsHaca = hacaSolver.solve()

print(f"Instanz {TAG}:")
print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG), gmmSolver.calculateAccurcacy(segmentsGmm, TAG), 
     acaSolver.calculateAccurcacy(segmentsAca, TAG), hacaSolver.calculateAccurcacy(segmentsHaca, TAG))
print(klClusterSolver.calculateMacroPrecision(segmentsKlCluster, TAG), gmmSolver.calculateMacroPrecision(segmentsGmm, TAG), 
     acaSolver.calculateMacroPrecision(segmentsAca, TAG), hacaSolver.calculateMacroPrecision(segmentsHaca, TAG))
print(klClusterSolver.calculateMacroRecall(segmentsKlCluster, TAG), gmmSolver.calculateMacroRecall(segmentsGmm, TAG), 
     acaSolver.calculateMacroRecall(segmentsAca, TAG), hacaSolver.calculateMacroRecall(segmentsHaca, TAG))
klClusterSolver.plotSegmentsAndGT([segmentsKlCluster, segmentsGmm, segmentsAca, segmentsHaca], ["KlCluster", "Gmm", "Aca", "Haca"])