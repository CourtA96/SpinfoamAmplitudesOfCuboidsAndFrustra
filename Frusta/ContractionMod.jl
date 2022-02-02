"""
julia> using ContractionMod
julia> j = [2, 2, 2];
julia> ints = ContractionMod.loadInts(j);
julia> sints = makeSizedArray.(ints);
julia> @time ContractionMod.Contraction(sints...)
  2.701165 seconds (11 allocations: 384 bytes)
7.900010707979186e-9 + 3.158254662302778e-24im

julia> using ContractionMod
julia> j = [3, 3, 3];
julia> ints = ContractionMod.loadInts(j);
julia> sints = makeSizedArray.(ints);
julia> @time ContractionMod.Contraction(sints...)
  6.323476 seconds (9 allocations: 160 bytes)
6.782085880953364e-11 - 8.310682282434477e-19im
"""
module ContractionMod

using LinearAlgebra
using StaticArrays

export fromTxt1 # Function to read in intertwiners
function fromTxt1(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[1],j[1],j[1],j[1],j[1],j[1]),[6,5,4,3,2,1])
    return arr
end


export fromTxt2 # Function to read in intertwiners
function fromTxt2(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[3],j[3],j[1],j[3],j[3],j[2]),[6,5,4,3,2,1])
    return arr
end


export fromTxt3 # Function to read in intertwiners
function fromTxt3(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[1],j[3],j[3],j[2],j[3],j[3]),[6,5,4,3,2,1])
    return arr
end


export fromTxt4 # Function to read in intertwiners
function fromTxt4(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[3],j[1],j[3],j[3],j[2],j[3]),[6,5,4,3,2,1])
    return arr
end


export fromTxt5 # Function to read in intertwiners
function fromTxt5(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[2],j[2],j[2],j[2],j[2],j[2]),[6,5,4,3,2,1])
    return arr
end


export fromTxt6 # Function to read in intertwiners
function fromTxt6(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[3],j[3],j[2],j[3],j[3],j[1]),[6,5,4,3,2,1])
    return arr
end


export fromTxt7 # Function to read in intertwiners
function fromTxt7(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[2],j[3],j[3],j[1],j[3],j[3]),[6,5,4,3,2,1])
    return arr
end


export fromTxt8 # Function to read in intertwiners
function fromTxt8(name::String, j::Vector{Int})

    io = open(name,"r")
    s = read(io,String)
    close(io)

    s = split(s)
    s = String.(s)
    sz = Int(size(s)[1]/3)
    s = reshape(s,3,sz)
    s = s[1,:].*s[2,:].*s[3,:]
    vec = [parse(Complex{Float64},ss) for ss in s]

    arr = permutedims(reshape(vec,j[3],j[2],j[3],j[3],j[1],j[3]),[6,5,4,3,2,1])
    return arr
end


export loadInts # function to load all intertwiners with given j at once
function loadInts(j::Vector) # j is a vector such as [2,2,2] or [3,3,3]

    #name = "Int1j" * string(j) * ".txt"
    name = "Int1j" * string("[",j[1],", ",j[1],", ",j[1],"]") * ".txt"
    #int1 = ContractionMod.fromTxt(name,j.+1)
    int1 = ContractionMod.fromTxt1(name,[j[1],j[1],j[1]].+1)

    #name = "Int2j" * string(j) * ".txt"
    #int2 = ContractionMod.fromTxt(name,j.+1)
    name = "Int2j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int2 = ContractionMod.fromTxt2(name,[j[1],j[2],j[3]].+1)

    #name = "Int3j" * string(j) * ".txt"
    #int3 = ContractionMod.fromTxt(name,j.+1)
    name = "Int3j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int3 = ContractionMod.fromTxt3(name,[j[1],j[2],j[3]].+1)

    name = "Int4j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int4 = ContractionMod.fromTxt4(name,[j[1],j[2],j[3]].+1)

    name = "Int5j" * string("[",j[2],", ",j[2],", ",j[2],"]") * ".txt"
    int5 = ContractionMod.fromTxt5(name,[j[2],j[2],j[2]].+1)

    name = "Int6j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int6 = ContractionMod.fromTxt6(name,[j[1],j[2],j[3]].+1)

    name = "Int7j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int7 = ContractionMod.fromTxt7(name,[j[1],j[2],j[3]].+1)

    name = "Int8j" * string("[",j[1],", ",j[2],", ",j[3],"]") * ".txt"
    int8 = ContractionMod.fromTxt8(name,[j[1],j[2],j[3]].+1)

    return int1,int2,int3,int4,int5,int6,int7,int8

end

export magIndices2
function magIndices2(dim::Int)
    J2 = dim-1
    return ntuple(i->J2-2*(i-1), dim)
end

export Contraction
function Contraction(int1::AbstractArray{Complex{Float64}, 6},
                     int2::AbstractArray{Complex{Float64}, 6},
                     int3::AbstractArray{Complex{Float64}, 6},
                     int4::AbstractArray{Complex{Float64}, 6},
                     int5::AbstractArray{Complex{Float64}, 6},
                     int6::AbstractArray{Complex{Float64}, 6},
                     int7::AbstractArray{Complex{Float64}, 6},
                     int8::AbstractArray{Complex{Float64}, 6})

    SUMre = Threads.Atomic{Float64}(0)
    SUMim = Threads.Atomic{Float64}(0)

    # Find the dimensions of each of the intertwiners

    dim1 = size(int1)
    dim2 = size(int2)
    dim3 = size(int3)
    dim4 = size(int4)
    dim5 = size(int5)
    dim6 = size(int6)
    dim7 = size(int7)

    # Loop over the given indices (all indices except con1, a6, con2, b4, con3, c3, and con4
    # which can be found in terms of the other indices).

    # @inbounds for a5 = 1:dim1[6], a4 = 1:dim1[5], a3 = 1:dim1[4], a2 = 1:dim1[3], a1 = 1:dim1[2]
    niters = dim1[6] * dim1[5] * dim1[4] * dim1[3] * dim1[2]
    @inbounds Threads.@threads for iter in 0:niters-1
        a5 = 1 + div(iter,(dim1[5] * dim1[4] * dim1[3] * dim1[2]))
        a4 = 1 + div(iter,(dim1[4] * dim1[3] * dim1[2])) % dim1[5]
        a3 = 1 + div(iter,(dim1[3] * dim1[2])) % dim1[4]
        a2 = 1 + div(iter,dim1[2]) % dim1[3]
        a1 = 1 + iter % dim1[2]

        SUM1 = convert(Complex{Float64},ContractionInt(int1, int2, int3, int4, int5, int6, int7, int8,
                              a1, a2, a3, a4, a5))

        Threads.atomic_add!(SUMre, real(SUM1))
        Threads.atomic_add!(SUMim, imag(SUM1))

    end

    SUM = Complex(SUMre[], SUMim[])
    return SUM
end

function ContractionInt(int1::AbstractArray{Complex{Float64}, 6},
                        int2::AbstractArray{Complex{Float64}, 6},
                        int3::AbstractArray{Complex{Float64}, 6},
                        int4::AbstractArray{Complex{Float64}, 6},
                        int5::AbstractArray{Complex{Float64}, 6},
                        int6::AbstractArray{Complex{Float64}, 6},
                        int7::AbstractArray{Complex{Float64}, 6},
                        int8::AbstractArray{Complex{Float64}, 6},
                        a1, a2, a3, a4, a5)

    SUM = zero(Complex{Float64})

    @inbounds begin

        # Find the dimensions of each of the intertwiners

        dim1 = size(int1)
        dim2 = size(int2)
        dim3 = size(int3)
        dim4 = size(int4)
        dim5 = size(int5)
        dim6 = size(int6)
        dim7 = size(int7)

        # find all of the magnetic indices

        #Mmg12 = magIndices2(dim1[2])
        #Mmg22 = magIndices2(dim1[3])
        #Mmg32 = magIndices2(dim1[1])
        #Mmg42 = magIndices2(dim3[2])
        #Mmg52 = magIndices2(dim2[2])
        #Mmg62 = magIndices2(dim3[1])

        Mmg1 = magIndices2(dim1[1])
        Mmg2 = magIndices2(dim5[1])
        Mmg3 = magIndices2(dim2[2])

        # Loop over the given indices (all indices except con1, a6, con2, b4, con3, c3, and con4
        # which can be found in terms of the other indices).

        # find magnetic index Con1
        Con12 = -Mmg1[a1] - Mmg1[a2] + Mmg1[a3] + Mmg1[a4] + Mmg1[a5]
        Jay2 = dim1[1]-1 # Value of j at index corresponging to A3
        # If the value of Con1 is allowed, convert to index notation and find A6
        if abs(Con12) <= Jay2 && (Con12 - Jay2) % 2 == 0
            y02 = 2+Jay2 # constant needed to convert to index notation
            con1 = div((-Con12+y02),2) # Index notation

            for a10 = 1:dim2[6], a9 = 1:dim2[5], a8 = 1:dim2[3], a7 = 1:dim2[2]

                # a6 = a7 - a8 + con1 - a9 + a10
                A62 = Mmg3[a7] - Mmg3[a8] + Mmg1[con1] - Mmg3[a9] + Mmg3[a10]
                Jay2 = dim2[1]-1
                # If the value of A6 is allowed, convert to index notation and find Con2
                if abs(A62) <= Jay2 && (A62 - Jay2) % 2 == 0
                    y02 = 2+Jay2
                    a6 = div((-A62+y02),2) # index notation

                    for b3 = 1:dim3[5], b2 = 1:dim3[1], b1 = 1:dim3[3]

                        # con2 = a8 - b1 - b2 + b3 + a2
                        Con22 = Mmg3[a8] - Mmg2[b1] - Mmg3[b2] + Mmg3[b3] + Mmg1[a2]
                        Jay2 = dim3[2]-1
                        # If the value of Con2 is allowed, convert to index notation and find B4
                        if abs(Con22) <= Jay2 && (Con22 - Jay2) % 2 == 0
                            y02 = 2+Jay2
                            con2 = div((-Con22+y02),2) # index notation

                            for b6 = 1:dim4[6], b5 = 1:dim4[4]

                                # b4 = -a7 - b6 + b5 + a1 - con2
                                B42 = -Mmg3[a7] - Mmg3[b6] + Mmg3[b5] + Mmg1[a1] + Mmg3[con2]
                                Jay2 = dim4[2]-1
                                # If the value of B4 is allowed, convert to index notation and find Con3
                                if abs(B42) <= Jay2 && (B42 - Jay2) % 2 == 0
                                    y02 = 2+Jay2
                                    b4 = div((-B42+y02),2) # index notation

                                    FACT1 = (int1[con1,a1,a2,a3,a4,a5] *
                                             int2[a6,a7,a8,con1,a9,a10] *
                                             int3[a8,con2,b1,b2,b3,a2] *
                                             int4[a7,b4,con2,b5,a1,b6])

                                    for c2 = 1:dim5[3], c1 = 1:dim5[2]

                                        # con3 = -c1 - c2 + a6 + b4 + b1
                                        Con32 = -Mmg2[c1] - Mmg2[c2] + Mmg2[a6] + Mmg2[b4] + Mmg2[b1]
                                        Jay2 = dim5[1]-1
                                        # If the value of Con3 is allowed, convert to index notation and find C4
                                        if abs(Con32) <= Jay2 && (Con32 - Jay2) % 2 == 0
                                            y02 = 2+Jay2
                                            con3 = div((-Con32+y02),2) # index notation

                                            FACT2 = FACT1 * int5[a6,c1,b1,con3,b4,c2]

                                            for c4 = 1:dim7[4]

                                                # c3 = a3 + b5 + c4 - con3 - b2
                                                C32 = Mmg1[a3] + Mmg3[b5] + Mmg3[c4] - Mmg2[con3] - Mmg3[b2]
                                                Jay2 = dim6[5]-1
                                                # If the value of C3 is allowed, convert to index notation and find Con4
                                                if abs(C32) <= Jay2 && (C32 - Jay2) % 2 == 0
                                                    y02 = 2+Jay2
                                                    c3 = div((-C32+y02),2) # index notation

                                                    # con4 = -a10 - a5 + c4 + b6 + c2
                                                    Con42 = -Mmg3[a10] - Mmg1[a5] + Mmg3[c4] + Mmg3[b6] + Mmg2[c2]
                                                    Jay2 = dim7[2]-1
                                                    # If the value of C3 is allowed, convert to index notation and contract
                                                    if abs(Con42) <= Jay2 && (Con42 - Jay2) % 2 == 0
                                                        y02 = 2+Jay2
                                                        con4 = div((-Con42+y02),2) # index notation

                                                        # perform the contraction by multiplying the specified entries of the intertwiner
                                                        # together and then adding the result to SUM.
                                                        SUM += (FACT2 *
                                                                int6[a3,b5,b2,con3,c3,c4] *
                                                                int7[a10,b6,a5,c4,con4,c2] *
                                                                int8[a9,a4,b3,c3,c1,con4])

                                                    end

                                                end

                                            end

                                        end

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    return SUM

end


function naiveContraction(int1::AbstractArray{Complex{Float64}, 6},
                        int2::AbstractArray{Complex{Float64}, 6},
                        int3::AbstractArray{Complex{Float64}, 6},
                        int4::AbstractArray{Complex{Float64}, 6},
                        int5::AbstractArray{Complex{Float64}, 6},
                        int6::AbstractArray{Complex{Float64}, 6},
                        int7::AbstractArray{Complex{Float64}, 6},
                        int8::AbstractArray{Complex{Float64}, 6})

    SUM = zero(Complex{Float64})

    # Find the dimensions of each of the intertwiners

    #int1 = big.(int1)
    #int2 = big.(int2)
    #int3 = big.(int3)
    #int4 = big.(int4)
    #nt5 = big.(int5)
    #int6 = big.(int6)
    #int7 = big.(int7)
    #int8 = big.(int8)

    dim1 = size(int1)
    dim2 = size(int2)
    dim3 = size(int3)
    dim4 = size(int4)
    dim5 = size(int5)
    dim6 = size(int6)
    dim7 = size(int7)

        # find all of the magnetic indices

        #Mmg12 = magIndices2(dim1[2])
        #Mmg22 = magIndices2(dim1[3])
        #Mmg32 = magIndices2(dim1[1])
        #Mmg42 = magIndices2(dim3[2])
        #Mmg52 = magIndices2(dim2[2])
        #Mmg62 = magIndices2(dim3[1])

        #Mmg1 = magIndices2(dim1[1])
        #Mmg2 = magIndices2(dim5[1])
        #Mmg3 = magIndices2(dim2[2])

        # Loop over the given indices (all indices except con1, a6, con2, b4, con3, c3, and con4
        # which can be found in terms of the other indices).

        # find magnetic index Con1
        #Con12 = -Mmg1[a1] - Mmg1[a2] + Mmg1[a3] + Mmg1[a4] + Mmg1[a5]
        #Jay2 = dim1[1]-1 # Value of j at index corresponging to A3
        # If the value of Con1 is allowed, convert to index notation and find A6
        #if abs(Con12) <= Jay2 && (Con12 - Jay2) % 2 == 0
        #    y02 = 2+Jay2 # constant needed to convert to index notation
        #    con1 = div((-Con12+y02),2) # Index notation

        for con1 = 1:dim1[1], a1 = 1:dim1[2], a2 = 1:dim1[3], a3 = 1:dim1[4], a4 = 1:dim1[5], a5 = 1:dim1[6],
            a6 = 1:dim2[1], a10 = 1:dim2[6], a9 = 1:dim2[5], a8 = 1:dim2[3], a7 = 1:dim2[2]

                # a6 = a7 - a8 + con1 - a9 + a10
                #A62 = Mmg3[a7] - Mmg3[a8] + Mmg1[con1] - Mmg3[a9] + Mmg3[a10]
                #Jay2 = dim2[1]-1
                # If the value of A6 is allowed, convert to index notation and find Con2
                #if abs(A62) <= Jay2 && (A62 - Jay2) % 2 == 0
                #    y02 = 2+Jay2
                #    a6 = div((-A62+y02),2) # index notation

            for con2 in 1:dim3[2], b3 = 1:dim3[5], b2 = 1:dim3[1], b1 = 1:dim3[3]

                        # con2 = a8 - b1 - b2 + b3 + a2
                        #Con22 = Mmg3[a8] - Mmg2[b1] - Mmg3[b2] + Mmg3[b3] + Mmg1[a2]
                        #Jay2 = dim3[2]-1
                        # If the value of Con2 is allowed, convert to index notation and find B4
                        #if abs(Con22) <= Jay2 && (Con22 - Jay2) % 2 == 0
                        #    y02 = 2+Jay2
                        #    con2 = div((-Con22+y02),2) # index notation

                for b4 in 1:dim4[2], b6 = 1:dim4[6], b5 = 1:dim4[4]

                                # b4 = -a7 - b6 + b5 + a1 - con2
                                #B42 = -Mmg3[a7] - Mmg3[b6] + Mmg3[b5] + Mmg1[a1] + Mmg3[con2]
                                #Jay2 = dim4[2]-1
                                # If the value of B4 is allowed, convert to index notation and find Con3
                                #if abs(B42) <= Jay2 && (B42 - Jay2) % 2 == 0
                                #    y02 = 2+Jay2
                                #    b4 = div((-B42+y02),2) # index notation

                            FACT1 = (int1[con1,a1,a2,a3,a4,a5] *
                                int2[a6,a7,a8,con1,a9,a10] *
                                int3[a8,con2,b1,b2,b3,a2] *
                                int4[a7,b4,con2,b5,a1,b6])

                            for con3 = 1:dim5[4], c2 = 1:dim5[3], c1 = 1:dim5[2]

                                        # con3 = -c1 - c2 + a6 + b4 + b1
                                        #Con32 = -Mmg2[c1] - Mmg2[c2] + Mmg2[a6] + Mmg2[b4] + Mmg2[b1]
                                        #Jay2 = dim5[1]-1
                                        # If the value of Con3 is allowed, convert to index notation and find C4
                                        #if abs(Con32) <= Jay2 && (Con32 - Jay2) % 2 == 0
                                        #    y02 = 2+Jay2
                                        #    con3 = div((-Con32+y02),2) # index notation

                                FACT2 = FACT1 * int5[a6,c1,b1,con3,b4,c2]

                                for c4 = 1:dim7[4], c3 = 1:dim6[5], con4 = 1:dim7[5]

                                                # c3 = a3 + b5 + c4 - con3 - b2
                                                #C32 = Mmg1[a3] + Mmg3[b5] + Mmg3[c4] - Mmg2[con3] - Mmg3[b2]
                                                #Jay2 = dim6[5]-1
                                                # If the value of C3 is allowed, convert to index notation and find Con4
                                                #if abs(C32) <= Jay2 && (C32 - Jay2) % 2 == 0
                                                #    y02 = 2+Jay2
                                                #    c3 = div((-C32+y02),2) # index notation

                                                    # con4 = -a10 - a5 + c4 + b6 + c2
                                                #    Con42 = -Mmg3[a10] - Mmg1[a5] + Mmg3[c4] + Mmg3[b6] + Mmg2[c2]
                                                #    Jay2 = dim7[2]-1
                                                #    # If the value of C3 is allowed, convert to index notation and contract
                                                #    if abs(Con42) <= Jay2 && (Con42 - Jay2) % 2 == 0
                                                #        y02 = 2+Jay2
                                                #        con4 = div((-Con42+y02),2) # index notation

                                                        # perform the contraction by multiplying the specified entries of the intertwiner
                                                        # together and then adding the result to SUM.

                                    SUM += (FACT2 *
                                        int6[a3,b5,b2,con3,c3,c4] *
                                        int7[a10,b6,a5,c4,con4,c2] *
                                        int8[a9,a4,b3,c3,c1,con4])

                    end

                end

            end

        end

    end

    return SUM

end


export makeSizedArray
function makeSizedArray(A::Array)
    SizedArray{Tuple{size(A)...}}(A)
end

function main2()
    j = [2, 2, 2]
    ints = ContractionMod.loadInts(j)
    sints = makeSizedArray.(ints)
    ContractionMod.Contraction(ints...)
end

function main3()
    j = [3, 3, 3]
    ints = ContractionMod.loadInts(j)
    sints = makeSizedArray.(ints)
    ContractionMod.Contraction(ints...)
end

end # module
