module createIntertwiners

using Cuba
using LinearAlgebra
using Distributed
using Memoize

export zr
function zr(num::Float64)
    if abs(num) < 2*(10^-7)
        return 0
    else
        return num
    end
end

@memoize function fact(j::Int,s::Int,M1::Float64,M2::Float64)
    return Float64(factorial(BigInt((j/2)+M1-s))*
            factorial(BigInt(s))*factorial(BigInt(M2-M1+s))*
            factorial(BigInt((j/2)-M2-s)))
end

@memoize function factCosSin(j::Int,s::Int,M1::Float64,M2::Float64,beta::Float64,dfact::Float64)
    return ((-1.0)^(s)*(1.0*im)^(M1-M2))*(dfact)^(-1)*
    cos(beta/2)^(j+M1-M2-2*s)*sin(beta/2)^(M2-M1+2*s)
end

export G
function G(alpha::Float64,beta::Float64,gamma::Float64,j::Int,M::Array{Any,2},d::Array{Complex{Float64},2})

    # Initialize matrices

    length = j+1 # size of matrix

    D = zeros(Complex{Float64},length,length)

    # Create Large D matrix from Wikipedia

    Dtemp = zeros(Complex{Float64},length,length)

    for I = 1:length, J= 1:length
        Dtemp[I,J] = cis(-M[I,J][2]*alpha-M[I,J][1]*gamma)
    end

    dee = zeros(Complex{Float64},length,length)
    dee[:,:] = d

    for I = 1:length, J = 1:length
       SUM = 0

        for s = 0:trunc(Int,j)
            if ((j/2)+M[I,J][1]-s)>=0 && (M[I,J][2]-M[I,J][1]+s)>=0 && ((j/2)-M[I,J][2]-s)>=0
                # If the value of s will not cause the factorials to be negative, find the
                # factorials and add the value to the sum

                dfact = fact(j,s,M[I,J][1],M[I,J][2])

                SUM += factCosSin(j,s,M[I,J][1],M[I,J][2],beta,dfact)

            end

        end

        # Multiply the sum by the square root

        dee[I,J] = dee[I,J]*SUM

        # Multiply the sum by the square root

    end

    # Combine the small d and large D matrices

    for I = 1:length, J = 1:length
        D[I,J] = Dtemp[I,J]*dee[I,J]
    end

    return D
end

# function to make a tensor with given parameters takes global rotation angles, and
# an array of strings that indicate direction and an array of j

export PreTensor
function PreTensor(j::Vector{Int},v::Vector{String},MM::Array{Any,1},DD::Array{Any,1})

    # Create a array that contains the length of each index

    L = Array{Int,1}(undef,3)

    for i = 1:3
        L[i] = trunc(Int,(j[i])+1)
    end

    # Find the longest length in L

    long = L[1]

    for i = 2:3
        if L[i] > long
            long = L[i]
        end
    end

    # Initialize and array that will store the parameters and an array that will hold
    # the initial z vector for each rotation

    A = zeros(Complex{Float64},long,3)

    z = zeros(Complex{Float64},long,3)

    for i = 1:3
        z[1,i] = 1.0 + 0*im
    end

    # Rotate the z vectors to match the input v and store in A.

    for i = 1:3
        if v[i] == "x"
            A[1:L[i],i] = G(-pi/2,pi/2,pi/2,j[i],MM[i],DD[i]) * z[1:L[i],i]
        elseif v[i] == "minusx"
            A[1:L[i],i] = G(-pi/2,-pi/2,pi/2,j[i],MM[i],DD[i]) * z[1:L[i],i]
        elseif v[i] == "y"
            A[1:L[i],i] = G(0.0,-pi/2,0.0,j[i],MM[i],DD[i]) * z[1:L[i],i]
        elseif v[i] =="minusy"
            A[1:L[i],i] = G(0.0,pi/2,0.0,j[i],MM[i],DD[i]) * z[1:L[i],i]
        elseif v[i] == "z"
            A[1:L[i],i] = z[1:L[i],i]
        elseif v[i] == "minusz"
            A[1:L[i],i] = G(0.0,pi*1.0,pi*1.0,j[i],MM[i],DD[i]) * z[1:L[i],i]
        else
            return error("InputError: Incorrect input, input must be x, y, z, minusx, minusy, or minusz")
        end
    end

    return A
end

export Tensor
function Tensor(alpha::Float64,beta::Float64,gamma::Float64,j::Vector{Int},
        MM::Array{Any,1},DD::Array{Any,1},A::Array{Complex{Float64},2})

        # Create a array that contains the length of each index

        L = Array{Int,1}(undef,3)

        for i = 1:3
            L[i] = trunc(Int,(j[i])+1)
        end

        # Create a table that will store the tensor values

        tbl = zeros(Complex{Float64},L[1],L[2],L[3],L[1],L[2],L[3]) # store the tensor as a 6D matrix

        # Calculate each of the wigner matrices

        G1 = G(alpha,beta,gamma,j[1],MM[1],DD[1])
        G1con = conj(transpose(G1))
        G2 = G(alpha,beta,gamma,j[2],MM[2],DD[2])
        G2con = conj(transpose(G2))
        G3 = G(alpha,beta,gamma,j[3],MM[3],DD[3])
        G3con = conj(transpose(G3))

        # pre-calculate the vectors

        v1 = (G1*A[1:L[1],1])
        v2 = (G2*A[1:L[2],2])
        v3 = (G3*A[1:L[3],3])
        v4 = ((A[1:L[1],1]')*G1con)
        v5 = (((A[1:L[2],2])')*G2con)
        v6 = ((A[1:L[3],3]')*G3con)

        for i = 1:L[1], J = 1:L[2], k = 1:L[3], l = 1:L[1], m = 1:L[2], n = 1:L[3]

            # Create the tensor by taking the product of g and a given column of A

            tbl[i,J,k,l,m,n] = v1[i]*v2[J]*v3[k]*v4[l]*v5[m]*v6[n]

        end

        return tbl

end

export kernel
function kernel(L1::Int,Inds1::Array{CartesianIndex{6},1},jay::Array{Int,1},MM::Array{Any,1},DD::Array{Any,1},A::Array{Complex{Float64},2})

    i = Inds1[L1][1]
    j = Inds1[L1][2]
    k = Inds1[L1][3]
    l = Inds1[L1][4]
    m = Inds1[L1][5]
    n = Inds1[L1][6]

    # integrate the integrand function and store in tbl

    result,err = cuhre((x, f) -> (f[1],f[2]) = reim((1/(16*(pi^2)))*sin(pi*x[2])*8*(pi^3)*
        Tensor((2*pi*x[1]),(pi*x[2]),(4*pi*x[3]),jay,MM,DD,A)[i,j,k,l,m,n]),3,2,minevals = 1e2,maxevals=1e4)
    return result[1]+result[2]*im

end

# Intertwiner() takes vectors as input and skips zeros and is parallelized
export Intertwiner
function Intertwiner(V::Array{String,1},Jay::Array{Int,1},MM::Array{Any,1},DD::Array{Any,1})
    # Create array of lengths of each index

    lngth = Array{Int,1}(undef,3)

    for i = 1:3
        lngth[i] = trunc(Int,(Jay[i])+1)
    end

    tbl1 = zeros(lngth[1],lngth[2],lngth[3],lngth[1],lngth[2],lngth[3])

    # store parameters as global variables so they can be called by the integrand function

    v = V
    jay = Jay

    A = PreTensor(jay,v,MM,DD)

    for i = 1:lngth[1], j = 1:lngth[2], k = 1:lngth[3], l = 1:lngth[1], m = 1:lngth[2]
        for n = 1:lngth[3]

            tbl1[i,j,k,l,m,n] = -(MM[1][i][1]+MM[2][j][1]+MM[3][k][1])+MM[1][l][1]+MM[2][m][1]+MM[3][n][1]

        end
    end

    Inds1 = findall(x->x==0,tbl1)

    tbl = zeros(Complex{Float64},lngth[1],lngth[2],lngth[3],lngth[1],lngth[2],lngth[3]) # Store the intertwiner in a 6D matrix

#        for L1 = 1:length(Inds1)
#            i = Inds1[L1][1]
#            j = Inds1[L1][2]
#            k = Inds1[L1][3]
#            l = Inds1[L1][4]
#            m = Inds1[L1][5]
#            n = Inds1[L1][6]
#            tbl1[i,j,k,l,m,n] = kernel(L1,Inds1,v,jay,MM,DD)
#        end

    ftbl = Array{Future}(undef, size(tbl))

    for L1 = 1:length(Inds1)
            i = Inds1[L1][1]
            j = Inds1[L1][2]
            k = Inds1[L1][3]
            l = Inds1[L1][4]
            m = Inds1[L1][5]
            n = Inds1[L1][6]
            p = workers()[mod1(L1, nworkers())]
            ftbl[i,j,k,l,m,n] = remotecall(kernel, p, L1,Inds1,jay,MM,DD,A)     # remember to use @everywhere before include
    end

    for L1 = 1:length(Inds1)
        i = Inds1[L1][1]
        j = Inds1[L1][2]
        k = Inds1[L1][3]
        l = Inds1[L1][4]
        m = Inds1[L1][5]
        n = Inds1[L1][6]
        tbl[i,j,k,l,m,n] = fetch(ftbl[i,j,k,l,m,n])
    end

    return tbl

end

export kernel2
function kernel2(L1::Int,Inds1::Array{CartesianIndex{6},1},jay::Array{Int,1},MM::Array{Any,1},DD::Array{Any,1},v::Array{String,1})

    i = Inds1[L1][1]
    j = Inds1[L1][2]
    k = Inds1[L1][3]
    l = Inds1[L1][4]
    m = Inds1[L1][5]
    n = Inds1[L1][6]

    # integrate the integrand function and store in tbl

    result,err = cuhre((x, f) -> (f[1],f[2]) = reim((1/(8*(pi^2)))*sin(pi*x[2])*4*(pi^3)*
        Tensor((2*pi*x[1]),(pi*x[2]),(2*pi*x[3]),jay,MM,DD,PreTensor(jay,v,MM,DD))[i,j,k,l,m,n]),3,2,minevals = 1e2,maxevals=1e4)
    return result[1]+result[2]*im

end

export Intertwiner2
function Intertwiner2(V::Array{String,1},Jay::Array{Int,1},MM::Array{Any,1},DD::Array{Any,1})
    # Create array of lengths of each index

    lngth = Array{Int,1}(undef,3)

    for i = 1:3
        lngth[i] = trunc(Int,(Jay[i])+1)
    end

    tbl1 = zeros(lngth[1],lngth[2],lngth[3],lngth[1],lngth[2],lngth[3])

    # store parameters as global variables so they can be called by the integrand function

    v = V
    jay = Jay

    for i = 1:lngth[1], j = 1:lngth[2], k = 1:lngth[3], l = 1:lngth[1], m = 1:lngth[2]
        for n = 1:lngth[3]

            tbl1[i,j,k,l,m,n] = -(MM[1][i][1]+MM[2][j][1]+MM[3][k][1])+MM[1][l][1]+MM[2][m][1]+MM[3][n][1]

        end
    end

    Inds1 = findall(x->x==0,tbl1)

    tbl = zeros(Complex{Float64},lngth[1],lngth[2],lngth[3],lngth[1],lngth[2],lngth[3]) # Store the intertwiner in a 6D matrix

    for L1 = 1:length(Inds1)
        i = Inds1[L1][1]
        j = Inds1[L1][2]
        k = Inds1[L1][3]
        l = Inds1[L1][4]
        m = Inds1[L1][5]
        n = Inds1[L1][6]
        tbl[i,j,k,l,m,n] = kernel2(L1,Inds1,jay,MM,DD,v)
    end

    return tbl

end

export generateMm
function generateMm(j::Vector{Int})

  L = Array{Int,1}(undef,3)

  for i = 1:3
      L[i] = trunc(Int,j[i]+1)
  end

  # L = ntuple(i -> j[i]+1, 3)

  long = j[1]+1
  for I = 2:3
      if j[I]+1 > long
          long = j[I]+1
      end
  end

  M1 = Array{Any,2}(undef,L[1],L[1])
  for I = 1:L[1], J = 1:L[1]
      M1[I,J] = [j[1]/2+(1-I),j[1]/2+(1-J)]
  end

  M2 = Array{Any,2}(undef,L[2],L[2])
  for I = 1:L[2], J = 1:L[2]
    M2[I,J] = [j[2]/2+(1-I),j[2]/2+(1-J)]
  end

  M3 = Array{Any,2}(undef,L[3],L[3])
  for I = 1:L[3], J = 1:L[3]
      M3[I,J] = [j[3]/2+(1-I),j[3]/2+(1-J)]
  end

  mm = Array{Any,1}(undef,3)
  mm[1] = M1
  mm[2] = M2
  mm[3] = M3

  D1 = zeros(Complex{Float64},L[1],L[1])

  for I = 1:L[1], J = 1:L[1]

      # First find the sqrt term outside the sum

      D1[I,J] = sqrt(factorial(trunc(BigInt,(j[1]/2)+M1[I,J][2]))*factorial(trunc(BigInt,(j[1]/2)-M1[I,J][2]))*
              factorial(trunc(BigInt,(j[1]/2)+M1[I,J][1]))*factorial(trunc(BigInt,(j[1]/2)-M1[I,J][1])))
  end

  D2 = zeros(Complex{Float64},L[2],L[2])

  for I = 1:L[2], J = 1:L[2]

      # First find the sqrt term outside the sum

      D2[I,J] = sqrt(factorial(trunc(BigInt,(j[2]/2)+M2[I,J][2]))*factorial(trunc(BigInt,(j[2]/2)-M2[I,J][2]))*
              factorial(trunc(BigInt,(j[2]/2)+M2[I,J][1]))*factorial(trunc(BigInt,(j[2]/2)-M2[I,J][1])))
  end

  D3 = zeros(Complex{Float64},L[3],L[3])

  for I = 1:L[3], J = 1:L[3]

      # First find the sqrt term outside the sum

      D3[I,J] = sqrt(factorial(trunc(BigInt,(j[3]/2)+M3[I,J][2]))*factorial(trunc(BigInt,(j[3]/2)-M3[I,J][2]))*
              factorial(trunc(BigInt,(j[3]/2)+M3[I,J][1]))*factorial(trunc(BigInt,(j[3]/2)-M3[I,J][1])))
  end

  dd = Array{Any,1}(undef,3)
  dd[1] = D1
  dd[2] = D2
  dd[3] = D3

  return [mm,dd]
end

export srt
function srt(arr::Array{Complex{Float64},6},iVec::Vector{Int})

    i1 = iVec[1]
    i2 = iVec[2]
    i3 = iVec[3]

    dim = size(arr)

    newArr = zeros(Complex{Float64}, dim[1],dim[2],dim[3],dim[4],dim[5],dim[6])

    for i = 1:dim[1], j = 1:dim[2], k = 1:dim[3], l = 1:dim[4], m = 1:dim[5], n = 1:dim[6]

        if i1==1 && i2==1 && i3==1
            newArr[i,j,k,l,m,n] = arr[i,j,k,l,m,n]
        elseif i1==2 && i2==1 && i3==1
            newArr[i,j,k,l,m,n] = arr[l,j,k,i,m,n]
        elseif i1==2 && i2==2 && i3==1
            newArr[i,j,k,l,m,n] = arr[l,m,k,i,j,n]
        elseif i1==2 && i2==2 && i3==2
            newArr[i,j,k,l,m,n] = arr[l,m,n,i,j,k]
        elseif i1==1 && i2==2 && i3==2
            newArr[i,j,k,l,m,n] = arr[i,j,k,l,m,n]
        elseif i1==1 && i2==1 && i3==2
            newArr[i,j,k,l,m,n] = arr[i,j,n,l,m,k]
        elseif i1==2 && i2==1 && i3==2
            newArr[i,j,k,l,m,n] = arr[l,j,n,i,m,k]
        elseif i1==1 && i2==2 && i3==1
            newArr[i,j,k,l,m,n] = arr[i,m,k,l,j,n]
        else
            return print("error")
        end

    end

    return newArr
end

export createInts
function createInts(Jay::Vector{Int})

    j1 = Jay[1]
    j2 = Jay[2]
    j3 = Jay[3]
    j4 = Jay[4]
    j5 = Jay[5]
    j6 = Jay[6]

    InterTemp = Array{Any}(undef,8,1)
    Inter = Array{Any}(undef,8,1)

    v = Array{Any}(undef,8,1)
        v[1] = ["x","y","z"]
        v[2] = ["x","minusy","z"]
        v[3] = ["minusx","y","z"]
        v[4] = ["x","y","minusz"]
        v[5] = ["minusx","y","minusz"]
        v[6] = ["x","y","minusz"]
        v[7] = ["x","minusy","z"]
        v[8] = ["minusx","y","z"]

    j = Array{Any}(undef,8,1)
        j[1] = [j3,j1,j2];
        j[2] = [j3,j5,j6];
        j[3] = [j6,j4,j2];
        j[4] = [j5,j1,j4];
        j[5] = [j3,j1,j2];
        j[6] = [j3,j5,j6];
        j[7] = [j6,j4,j2];
        j[8] = [j5,j1,j4];


    srtOrder = Array{Any}(undef,8,1)
        srtOrder[1] = [1,1,1]
        srtOrder[2] = [1,2,1]
        srtOrder[3] = [2,1,1]
        srtOrder[4] = [1,1,2]
        srtOrder[5] = [2,1,2]
        srtOrder[6] = [1,1,2]
        srtOrder[7] = [1,2,1]
        srtOrder[8] = [2,1,1]

        count = 0

    for i = 1:8

        name1 = "Int"
        name2 = string(i)
        name3 = "j"
        name4 = string(j[i])
        name5 = ".txt"
        name = name1*name2*name3*name4*name5

        if isfile(name) == false

            A = generateMm(j[i])
            Mm = A[1]
            Dd = A[2]

            print("\n The time to compute intertwiner ",string(i)," for j = ",string(j[i])," is: ")

            @time InterTemp[i] = Intertwiner(v[i],j[i],Mm,Dd)

            print("\n The time to sort the intertwiner is: ")

            @time Inter[i] = srt(InterTemp[i],srtOrder[i])

            toTxt(Inter[i],name)

            print("\n The norm of the intertwiner is: ",
            sum(abs.(Inter[i])))
            print("\n\n")

            count += 1

        end

    end

    print("\n The number of intertwiners computed is: ", string(count),"\n\n")


end

export toTxt
function toTxt(inter::Array{Complex{Float64},6},name::String)

    dims = size(inter)
    inter = string.(inter)

    io = open(name,"w")
    for i = 1:dims[1],j = 1:dims[2],k = 1:dims[3],l=1:dims[4],m=1:dims[5],n=1:dims[6]
        write(io,inter[i,j,k,l,m,n])
        write(io, " ")
    end
    close(io)
end

export test
function test(jay::Array{Int,1})

    v = ["x","y","z"]
    A = generateMm(jay)
    MM = A[1]
    DD = A[2]

    A = PreTensor(jay,v,MM,DD)

    integ = cuhre((x, f) -> (f[1],f[2]) = reim((1/(8*(pi^2)))*sin(pi*x[2])*4*(pi^3)*
        Tensor((2*pi*x[1]),(pi*x[2]),(2*pi*x[3]),jay,MM,DD,A)[1,1,1,1,1,1]),3,2,minevals = 100,maxevals=1e4)

    return integ

end

export test2
function test2()

    jay = [1,1,1]
    v = ["x","y","z"]
    A = generateMm(jay)
    MM = A[1]
    DD = A[2]

    integ = Intertwiner(v,jay,MM,DD)

    return integ

end

export test3
function test3()

    jay = [1,1,1]
    v = ["x","y","z"]
    A = generateMm(jay)
    MM = A[1]
    DD = A[2]

    integ = Intertwiner2(v,jay,MM,DD)

    return integ

end

end # module
