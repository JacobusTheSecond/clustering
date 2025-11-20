import os

from drifterSolvers import KlClusterDriftersSolver


klcluster = True
datafolder = "/home/charon/data/drifter"
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]
lowerbounds = []
solutionsizes = []


solver = KlClusterDriftersSolver(drifterFiles[0:10],freedelta=4,simpdelta=0.5, withaggressive=False,rossbyRadius=False) if klcluster else None
    #lowerbounds.append(solver.lowerbound(False))
#print(f"Size of Sol: {solver.solutionSize(False)}")

solver.solutionSize(False)

def printCinterval(i):
    print(i.curve,i.start.value,i.end.value,"(",solver.cc.mapToBase(i.curve,i.start).value,solver.cc.mapToBase(i.curve,i.end).value,")")
for c in solver.clusters:
    printCinterval(c.center())
    for m in c:
        print(" ",end="")
        printCinterval(m)
#solver.drawInputField()
#solver.drawResultField()
#solver.plotInputAndResult()
#solver.plotCurve()
#solver.plotResults()
