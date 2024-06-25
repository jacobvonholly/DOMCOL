using JuMP
using HiGHS
using SparseArrays
using Timers
model = read_from_file("supportcase12.mps")
data = lp_matrix_data(model);


global dom_pairs = []
global A = data.A
global c = data.c
global n = length(c)
global m = floor(Int, length(A) / n)
global dom_pairs = []
global zbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


function CreateBitstring(x)
    local varbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if issparse(A[:, x])
        for i in findnz(A[:, x])[1]
            if i % 64 == 0
                varbitstring[64] = 1
            else
                varbitstring[i % 64] = 1
            end
        end
        return(varbitstring)
    end
    for i in 1 : 64
        local p = 0
        while (p*64 + i) <= m && A[(p*64 + i), x] == 0
            p += 1
        end
        if (p*64 + i) <= m
            varbitstring[i] = 1
        end
    end
    return(varbitstring)
end


function findHeap(bitss)
    l = length(L)
    upper = l + 1
    lower = 0
    while upper - lower > 1
        pos = floor(Int, (upper + lower)/2)
        if L[pos][1] == bitss
            return(true, pos)
        elseif L[pos][1] < bitss
            lower = pos
        elseif L[pos][1] > bitss
            upper = pos
        end
    end
    return(false, upper)
end


function CreateHeaps()
    global L = []
    for i in 1 : n
        temp = CreateBitstring(i)
        # println(i)
        # println(i, " with bitstring ", temp)
        (inside, loc) = findHeap(temp)
        if inside
            push!(L[loc][2], i)
            # println(" now in same Heap with bitstring ", L[loc][1])
        else
            insert!(L, loc, [temp,[i]])
        end
    end
    return(L)
end


function findpairs(L)
    local l = length(L)
    local i1 = 1
    while i1 < l
        local i2 = i1 + 1
        while i2 <= l
            compare(L[i1],L[i2])
            i2 += 1
        end
        # print("nothing on ", string(i1))
        i1 += 1
    end
end

function compare(i1, i2)
    local p = 1
    local equal = 0
    local comp = ==
    if c[i1] < c[i2]
        comp = <=
    end
    if c[i1] > c[i2]
        comp = >=
    end
    if issparse(A[:, i1]) && issparse(A[:, i2])
        indices = sort(unique(append!(findnz(A[:, i1])[1], findnz(A[:, i2])[1])))
        leni = length(indices)
        while p <= leni && A[indices[p], i1] == A[indices[p], i2] && comp == ==
            p += 1
        end
        if p <= leni && A[indices[p], i1] < A[indices[p], i2]
            if comp == >=
                return()
            end
            comp = <=
        end
        if p <= leni && A[indices[p], i1] > A[indices[p], i2]
            if comp == <=
                return()
            end
            comp = >=
        end
        while p <= leni && comp(A[indices[p], i1], A[indices[p], i2])
            p += 1
        end
        if p == leni + 1
            push!(dom_pairs, [i1, i2])
            # print(i1, i2)
        end
    else
        while A[p, i1] == A[p, i2] && p <= m && comp == ==
            p += 1
        end
        if A[p, i1] < A[p, i2]
            if comp == >=
                return()
            end
            comp = <=
        end
        if A[p, i1] > A[p, i2]
            if comp == <=
                return()
            end
            comp = >=
        end
        while p <= m && comp(A[p, i1], A[p, i2])
            p += 1
        end
        if p == m + 1
            push!(dom_pairs, [i1, i2])
            # println(i1, i2)
        end
    end
end

#global L
tic()
CreateHeaps()
println("Creating Heaps took: ")
toc()
tic()
test = []
#=for i in 1 : n
    append!(test, i)
end 
L = [[1,test]] =#
l = length(L)
for i in 1 : l
    i1 = 1
    templ = length(L[i][2])
    while i1 < templ
        i2 = i1 + 1
        while i2 <= templ
            # println("comparing ", L[i][2][i1], " with ", L[i][2][i2])
            compare(L[i][2][i1], L[i][2][i2])
            i2 += 1
        end
        # println("that was all for ", L[i][2][i1])
        i1 += 1
    end
    #println(i/l, "%")
end
println("Runtime was: ")
toc()
println(" there are this many pairs: ", length(dom_pairs))
            










