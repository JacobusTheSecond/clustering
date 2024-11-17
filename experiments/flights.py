from flightSolvers import KlClusterFlightSolver

n = 300

solver = KlClusterFlightSolver(n)
solver.solve()

#solver.plotInput()
solver.plotInputAndResult()
solver.plotResults()