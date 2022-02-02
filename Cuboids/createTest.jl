using Distributed#, ClusterManagers
#addprocs(32)
#addprocs_slurm(64)
@everywhere include("createIntertwiners.jl")

#j = Array{Any}(undef,5,1)
#    j[1] = [1,1,1,1,1,1]
#    j[2] = [2,1,1,1,1,1]
#    j[3] = [3,1,1,1,1,1]
#    j[4] = [4,1,1,1,1,1]
#    j[5] = [5,1,1,1,1,1]

print("\nThe number of workers is: ",string(nworkers()),"\n")

createIntertwiners.createInts([1,1,1,1,1,1])

#for i = 1:5

#    print("\n********** New j **********\n")

#    createIntertwiners.createInts(j[i]) # createInts computes all of the eight intertwiners for the given j
                                        # it also checks to see if an intertwiner already exists in the file
                                        # so it only computes the necessary intertwiners

#end
