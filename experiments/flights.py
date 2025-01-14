from flightSolvers import KlClusterFlightSolver

#n = 300
ns = [25,50,100,200,400,800]
lowerbounds = []
solutionsizes = []

for n in ns:
    solver = KlClusterFlightSolver(n)
    lowerbounds.append(solver.lowerbound(False))
    solutionsizes.append(solver.solutionSize(False))
    print(f"ns: {ns}")
    print(f"Lowerbounds: {lowerbounds}")
    print(f"Size of Sol: {solutionsizes}")

#solver.plotInput()
solver.plotInputAndResult()
solver.plotResults()