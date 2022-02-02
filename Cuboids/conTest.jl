include("ContractionMod.jl")

Jvector = [1,1,1,1,1,1]

int1,int2,int3,int4,int5,int6,int7,int8 = ContractionMod.loadInts(Jvector)

println()

print("\n The time to compute the contraction for ",string(Jvector)," is: ")

@time contract = ContractionMod.Contraction(int1,int2,int3,int4,int5,int6,int7,int8)

print("\n The value of the contraction for ",string(Jvector)," is: ")

display(contract)

println()
