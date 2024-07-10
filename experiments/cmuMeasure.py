from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver
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
    tmmSolver = TmmCMUSolver(TAG)


    segmentsKlCluster = klClusterSolver.solve()
    #print(klClusterSolver.calculateAccurcacy(segmentsKlCluster, TAG))

    segmentsGmm = gmmSolver.solve()

    segmentsAca = acaSolver.solve()
    segmentsHaca = hacaSolver.solve()
    segmentsTmm = tmmSolver.solve()
    acc = {
        "KlCluster" : klClusterSolver.calculateAccuracy(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateAccuracy(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateAccuracy(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateAccuracy(segmentsHaca, TAG),
        "tmm" : klClusterSolver.calculateAccuracy(segmentsTmm, TAG),
    }
    trueAcc = {
        "KlCluster" : klClusterSolver.calculateAccuracyTrueLabels(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateAccuracyTrueLabels(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateAccuracyTrueLabels(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateAccuracyTrueLabels(segmentsHaca, TAG),
        "tmm" : klClusterSolver.calculateAccuracyTrueLabels(segmentsTmm, TAG),
    }
    macroPrec = {
        "KlCluster" : klClusterSolver.calculateMacroPrecision(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateMacroPrecision(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateMacroPrecision(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateMacroPrecision(segmentsHaca, TAG),
        "tmm" : klClusterSolver.calculateMacroPrecision(segmentsTmm, TAG),
    }
    macroRec = {
        "KlCluster" : klClusterSolver.calculateMacroRecall(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateMacroRecall(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateMacroRecall(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateMacroRecall(segmentsHaca, TAG),
        "tmm" : klClusterSolver.calculateMacroRecall(segmentsTmm, TAG),
    }
    macroF1 = {
        "KlCluster" : klClusterSolver.calculateMacroF1(segmentsKlCluster, TAG),
        "gmm" : klClusterSolver.calculateMacroF1(segmentsGmm, TAG),
        "aca" : klClusterSolver.calculateMacroF1(segmentsAca, TAG),
        "haca" : klClusterSolver.calculateMacroF1(segmentsHaca, TAG),
        "tmm" : klClusterSolver.calculateMacroF1(segmentsTmm, TAG),
    }
    result[TAG] = {"acc" : acc, "trueAcc" : trueAcc, "macroPrec" : macroPrec, "macroRec" : macroRec, "macroF1" : macroF1}
with open(file_path, 'w') as file:
    json.dump(result, file, indent=4)