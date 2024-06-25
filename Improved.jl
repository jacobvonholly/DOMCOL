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
global randvec = []

function genranvec()
    i = 1
    while i <= m+1
        temp = rand()
        push!(randvec, temp)
        i+=1
    end
end

genranvec()

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

function getIndices(bitss)
    I = []
    i = 1
    while i <= 64 
        if bitss[i] == 1
            j = 0
            while j*64 + i <= m
                append!(I, j*64 + i)
                j +=1
            end
        end
    end
    return(I)
end

function SparseIndices(X)
    I = []
    i = 1
    lenn = length(X[2])
    while i <= lenn
        append!(I, findnz(A[:, X[2][i]])[1])
        i += 1
    end
    return(unique(I))
end

function getrandaverage(x)
    a = randvec[m+1] * c[x]
    templist = findnz(A[:,x])
    lennn = length(templist[1])
    i = 1
    while i <= lennn
        a += randvec[templist[1][i]] * A[templist[1][i], x]
        i += 1
    end
    return(a)
end

function firstSort(List)
    # I = getIndices(List[1])
    # I = SparseIndices(List[2])
    len = length(List[2])
    templist = []
    i = 1
    while i <= len
        push!(templist, [Int(List[2][i]), [getrandaverage(List[2][i])]])
        i += 1
    end
    sort!(templist, by = x -> x[2][1])
    return(templist)
end

function dompairs(List)
    len = length(List[2])
    # I = getIndices(List[1])
    I = SparseIndices(List)
    global order = firstSort(List)
    global potdoms = []
    i = 1
    while i <= len
        temp = []
        j = 1
        while j < i
            push!(temp, j)
            j += 1
        end
        i += 1
        push!(potdoms, temp)
    end
    for i in I
        templistint = []
        j = 1
        while j <= len
            append!(templistint, j)
            j+= 1
        end
        dommergesort(templistint, i)
    end
    i = 1
    while i <= len
        lennn = length(potdoms[i])
        j = 1
        while j <= lennn
            push!(dom_pairs, [order[i][1], order[potdoms[i][j]][1]])
            # println(dom_pairs)
            j += 1
        end
        i += 1
    end
end

function dommergesort(ordering, y)
    # println("Calling dommergesort on: ", ordering, " and ", y)
    le = length(ordering)
    if le > 1
        return(dommerge(dommergesort(ordering[1:floor(Int, le/2)], y), dommergesort(ordering[floor(Int, le/2)+1:le], y),y))
    else
        return(ordering)
    end
end

function dommerge(order1, order2, y)
    # println("Calling dommerge on: ", order1, ", ", order2, " and ", y )
    neworder = []
    le1 = length(order1)
    le2 = length(order2)
    i1 = 1
    i2 = 1
    while i1 <= le1 && i2 <= le2
        if A[y, order[order1[i1]][1]] <= A[y, order[order2[i2]][1]]
            append!(neworder, order1[i1])
            i1 += 1
        else
            # println("order1[i1:le1]: ", order1[i1:le1])
            # println("potdoms[order2[i2]]: ", potdoms[order2[i2]])
            # println("filtered out: ", filter!(e-> (e in order1[i1:le1]), potdoms[order2[i2]]))
            filter!(e-> !(e in order1[i1:le1]), potdoms[order2[i2]])
            append!(neworder, order2[i2])
            i2 += 1
        end
    end
    if i1 > le1 
        append!(neworder, order2[i2:le2])
    else 
        append!(neworder, order1[i1:le1])
    end
    return(neworder)
end
tic()
CreateHeaps()
println("Creating Heaps took: ")
toc()
function main()
    i = 1
    lenngth = length(L)
    while i <= lenngth
        dompairs(L[i])
        i += 1
    end
end
tic()
main()
println("Runtime was: ")
toc()

println(" Anzahl: ", length(dom_pairs))
