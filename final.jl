using JuMP
using HiGHS
using SparseArrays
using Timers



function main()
    global runtime = 0
    global runtime2 = 0
    global start = time()
    global Averagevec = []
    global i = 1
    while i <= n
        append!(Averagevec, getrandaverage(i))
        global i += 1
    end
    i = 1
    lenngth = length(L)
    while i <= lenngth && runtime2 <= 1000
        global start2 = time()
        dompairs(L[i])
        i += 1
        global runtime2 += time() - start2
    end
    unique!(dom_pairs)
    global runtime = time() - start
    return(runtime,length(dom_pairs))
end

function Averagevect(i)
    if i > 0
        return(Averagevec[i])
    else
        return(-Averagevec[-i])
    end
end

function Av(i,j)
    if j > 0
        return(A[i, j])
    else
        return(-A[i, -j])
    end
end

function cv(i)
    if i > 0
        return(c[i])
    else
        return(-c[-i])
    end    
end


function genranvec()
    i = 1
    while i <= m+1
        temp = rand()
        push!(randvec, temp)
        i+=1
    end
end


function BetterHeaps()
    global L = []
    global LengthL = 0
    templist = []
    local templength = n
    genranvec()
    i = 1
    while i <= n
        push!(templist, [0, [i]])
        i += 1
    end
    sort!(templist)
    while templist != []
        genranvec()
        i = 1
        while i <= length(templist)
            templist[i][1] = getrandaverage(templist[i][2][1])
            i+=1
        end
        push!(L, templist[1:floor(Int, sqrt(templength))])
        LengthL += 1
        append!(L[LengthL], templist[max(floor(Int, sqrt(templength))+1,templength-floor(Int, sqrt(templength))):templength])
        templist = templist[floor(Int, sqrt(templength))+1:max(floor(Int, sqrt(templength))+1,templength-floor(Int, sqrt(templength)))-1]
        templength = length(templist)
    end
    i = 1
    while i <= LengthL
        templength = length(L[i])
        templist = []
        j = 1
        while j <= templength
            push!(templist, L[i][j][2][1])
            j += 1
        end
        L[i] = [0.0, templist]
        i += 1
    end
end

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
        push!(templist, [- Int(List[2][i]), [-Averagevec[List[2][i]]]])
        i += 1
    end
    sort!(templist, by = x -> x[2][1])
    return(templist)
end

function sortrows(Rows)
    olen = length(order)
    rlen = length(Rows)
    global Indorder = []
    i = 1
    while i <= rlen
        push!(Indorder, [rowsort(Rows[i]), [Rows[i]]])
        i += 1
    end    
    sort!(Indorder)
end

function rowsort(x)
    temprow = []
    olen = length(order)
    i = 1
    while i <= olen
        push!(temprow, [A[x, order[i][1]], i])
        i += 1
    end
    sort!(temprow)
    res = 0
    i = 1
    while i <= olen
        res += i*temprow[i][2]
        i += 1
    end
    return(res)
end

function dompairs(List)
    # I = getIndices(List[1])
    I = SparseIndices(List)
    global order = firstSort(List)
    len = length(order)
    # sortrows(I)
    global potdoms = []
    i = 1
    while i <= len
        temp = []
        j = 1
        while j < i
            if cv(order[j][1]) <= cv(order[i][1]) && ( order[j][1] > 0 || ( abs(order[j][1]) <  order[i][1])) && order[j][1] != -order[i][1]
                push!(temp, j)
            end
            j += 1
        end
        i += 1
        push!(potdoms, temp)
    end
    templistint = []
    j = 1
    while j <= len
        append!(templistint, j)
        j+= 1
    end
    i = 1
    lenr = length(I)
    global runtime3 = 0
    while i <= lenr && runtime3 <= 1000
        global start3 = time()
        # dommergesort(templistint, Indorder[i][2][1])
        dommergesort(templistint, I[i])
        global runtime3 += time() - start3
        i += 1
    end
    if i <= lenr
        return()
    end
    i = 1
    while i <= len
        lennn = length(potdoms[i])
        j = 1
        while j <= lennn
            push!(dom_pairs, [sign(order[i][1]order[potdoms[i][j]][1])min(abs(order[i][1]), abs(order[potdoms[i][j]][1])), max(abs(order[i][1]), abs(order[potdoms[i][j]][1]))])
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


function binary_search(array, target::Int64, starts::Int64)
    l = starts
    r = length(array)
    while l <= r
        pos = floor(Int, (l+r)/2)
        if array[pos] < target
            l = pos + 1
        elseif array[pos] > target
            r = pos - 1
        else
            return(true, pos)
        end
    end
    return(false, l)
end

function fastcut(L1, element)
    sort!(L1)
    l1 = length(L1)
    l2 = length(potdoms[element])
    if l1 == 0 || l2 == 0
        return()
    end
    if l1 < l2 
        i1 = 1
        u2 = 1
        while i1 <= l1 && u2 <= l2
            res = binary_search(potdoms[element], L1[i1], u2)
            if res[1]
                deleteat!(potdoms[element], res[2])
            end
            i1 += 1
            u2 = res[2]
        end
    else
        i2 = 1
        u1 = 1
        while i2 <= l2 && u1 <= l1
            res = binary_search(L1, potdoms[element][i2], u1)
            if res[1]
                deleteat!(potdoms[element], i2)
                i2 -= 1
                l2 -= 1
            end
            i2 += 1
            u1 = res[2]
        end
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
        if Av(y, order[order1[i1]][1]) <= Av(y, order[order2[i2]][1])
            append!(neworder, order1[i1])
            i1 += 1
        else
            # println("order1[i1:le1]: ", order1[i1:le1])
            # println("potdoms[order2[i2]]: ", potdoms[order2[i2]])
            # println("filtered out: ", filter!(e-> (e in order1[i1:le1]), potdoms[order2[i2]]))
            # filter!(e-> !(e in order1[i1:le1]), potdoms[order2[i2]])
            fastcut(order1[i1:le1], order2[i2])
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

function compareneg(i1, i2)
    local p = 1
    local equal = 0
    local comp = ==
    if -c[i1] < c[i2]
        comp = <=
    end
    if -c[i1] > c[i2]
        comp = >=
    end
    if issparse(A[:, i1]) && issparse(A[:, i2])
        indices = sort(unique(append!(findnz(A[:, i1])[1], findnz(A[:, i2])[1])))
        leni = length(indices)
        while p <= leni && -A[indices[p], i1] == A[indices[p], i2] && comp == ==
            p += 1
        end
        if p <= leni && -A[indices[p], i1] < A[indices[p], i2]
            if comp == >=
                return()
            end
            comp = <=
        end
        if p <= leni && -A[indices[p], i1] > A[indices[p], i2]
            if comp == <=
                return()
            end
            comp = >=
        end
        while p <= leni && comp(-A[indices[p], i1], A[indices[p], i2])
            p += 1
        end
        if p == leni + 1
            push!(dom_pairs, [-i1, i2])
            # print(i1, i2)
        end
    else
        while -A[p, i1] == A[p, i2] && p <= m && comp == ==
            p += 1
        end
        if -A[p, i1] < A[p, i2]
            if comp == >=
                return()
            end
            comp = <=
        end
        if -A[p, i1] > A[p, i2]
            if comp == <=
                return()
            end
            comp = >=
        end
        while p <= m && comp(-A[p, i1], A[p, i2])
            p += 1
        end
        if p == m + 1
            push!(dom_pairs, [-i1, i2])
            # println(i1, i2)
        end
    end
end

#test = []
#=for i in 1 : n
    append!(test, i)
end 
L = [[1,test]] =#
function Naivemain()
    global runtime = 0
    global start = time()
    global runtime2 = 0
    l = length(L)
    i = 1
    while i <= l && runtime2 <= 1000
        global start2 = time()
        i1 = 1
        templ = length(L[i][2])
        while i1 < templ
            i2 = i1 + 1
            while i2 <= templ
                # println("comparing ", L[i][2][i1], " with ", L[i][2][i2])
                compare(L[i][2][i1], L[i][2][i2])
                compareneg(L[i][2][i1], L[i][2][i2])
                i2 += 1
            end
            # println("that was all for ", L[i][2][i1])
            i1 += 1
        end
        #println(i/l, "%")
        i += 1
        global runtime2 += time() - start2
    end
    global runtime += time() - start
    return(runtime, length(dom_pairs))
end


cd("C:/Users/jacob/Desktop/Testing")

directory = readdir()

global directoryindex = 220
global const directorylength = length(directory)

while directoryindex <= directorylength
    println(time())
    global start = time()
    println(directory[directoryindex])
    model = read_from_file(directory[directoryindex])
    data = lp_matrix_data(model);

    global dom_pairs = []
    global A = data.A
    global c = data.c
    global n = length(c)
    global m = floor(Int, length(A) / n)
    global dom_pairs = []
    global zbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    global randvec = []
    genranvec()
    CreateHeaps()

    global preprocessingtime = time() - start

    result = main()
    println("one")
    improvedtime = result[1]
    improvedpairs = result[2]

    dom_pairs = []
    result = Naivemain()
    normaltime = result[1]
    normalpairs = result[2]

    a = 1000
    open("0000results.txt","a") do io
        println(io, "#", directory[directoryindex], "#",improvedpairs-normalpairs, "#", improvedtime , "#",  normaltime, "#", n, "#", m, "#", length(findnz(A)[3]), "#", improvedpairs, "#" )
    end

    println(directoryindex)

    global directoryindex += 1

end

global directoryindex = 3
global const directorylength = length(directory)

while directoryindex <= directorylength

    global start = time()

    model = read_from_file(directory[directoryindex])
    data = lp_matrix_data(model);

    global dom_pairs = []
    global A = data.A
    global c = data.c
    global n = length(c)
    global m = floor(Int, length(A) / n)
    global dom_pairs = []
    global zbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    global randvec = []
    genranvec()
    L = [1,[]]
    global i = 1
    while i <= n
        push!(L[2], i)
        global i += 1
    end

    preprocessingtime = time() - start

    result = main()
    improvedtime = result[1]
    improvedpairs = result[2]

    dom_pairs = []
    result = Naivemain()
    normaltime = result[1]
    normalpairs = result[2]

    a = 1000
    open("000fullresults.txt","a") do io
        println(io, "#", directory[directoryindex], "#",improvedpairs-normalpairs, "#", improvedtime , "#",  normaltime, "#", n, "#", m, "#", length(findnz(A)[3]), "#", improvedpairs, "#" )
    end

    println(directoryindex)

    global directoryindex += 1

end