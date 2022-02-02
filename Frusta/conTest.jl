include("ContractionMod.jl")

#j = Array{Any}(undef,5,1)
#    j[1] = [1,1,1,1,1,1]
#    j[2] = [2,1,1,1,1,1]
#    j[3] = [3,1,1,1,1,1]
#    j[4] = [4,1,1,1,1,1]
#    j[5] = [5,1,1,1,1,1]

Jvector = [1,1,1]

int1,int2,int3,int4,int5,int6,int7,int8 = ContractionMod.loadInts(Jvector)

println()

print("\n The time to compute the contraction for ",string(Jvector)," is: ")

@time contract = ContractionMod.Contraction(int1,int2,int3,int4,int5,int6,int7,int8)

print("\n The value of the contraction for ",string(Jvector)," is: ")

display(contract)

println()

#for i = 1:5

#    performContraction.contractInts(j[i]) # contractInts() contracts the eight intertwiners corresponding
                                          # to the given j

#end
