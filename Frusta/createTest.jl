using Distributed#, ClusterManagers

@everywhere include("createIntertwiners.jl")

print("\nThe number of workers is: ",string(nworkers()),"\n")

createIntertwiners.createInts([1,1,1])
