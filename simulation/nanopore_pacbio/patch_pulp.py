import pulp
if not hasattr(pulp, "list_solvers"):
    pulp.list_solvers = pulp.listSolvers
