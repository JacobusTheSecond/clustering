from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
import os
import json

file_path = "result.json"

result = {}

ITERATIONS = 10

if os.path.exists(file_path):
    with open(file_path, "r") as file:
        result = json.load(file)

for x in range(10, 15):
    TAG = x
    
    print(TAG, "TAG")
    klClusterSolver = KlClusterCMUSolver(TAG)
    gmmSolver = AcaCMUSolver(TAG, "gmm")
    acaSolver = AcaCMUSolver(TAG, "aca")
    hacaSolver = AcaCMUSolver(TAG, "haca")


    segmentsKlCluster = klClusterSolver.solve()
    #print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG))

    segmentsGmm = gmmSolver.solve()

    segmentsAca = acaSolver.solve()
    segmentsHaca = hacaSolver.solve()
    acc = {
        "KlCluster" : klClusterSolver.calculateAccuracy(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateAccuracy(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateAccuracy(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateAccuracy(segmentsHaca, TAG),
    }
    macroPrec = {
        "KlCluster" : klClusterSolver.calculateMacroPrecision(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateMacroPrecision(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateMacroPrecision(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateMacroPrecision(segmentsHaca, TAG),
    }
    macroRec = {
        "KlCluster" : klClusterSolver.calculateMacroRecall(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateMacroRecall(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateMacroRecall(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateMacroRecall(segmentsHaca, TAG),
    }
    result[TAG] = {"acc" : acc, "macroPrec" : macroPrec, "macroRec" : macroRec}
with open(file_path, 'w') as file:
    json.dump(result, file, indent=4)
