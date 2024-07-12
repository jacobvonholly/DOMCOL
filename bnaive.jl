using JuMP
using HiGHS
using SparseArrays
using Timers


function CreateBitstring(x)
    local varbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if issparse(A[:, x])
        for i in findnz(A[:, x])[1]
            if i % 32 == 0
                varbitstring[32] = 1
            else
                varbitstring[i % 32] = 1
            end
        end
        return(varbitstring)
    end
    for i in 1 : 32
        local p = 0
        while (p*32 + i) <= m && A[(p*32 + i), x] == 0
            p += 1
        end
        if (p*32 + i) <= m
            varbitstring[i] = 1
        end
    end
    return(varbitstring)
end

function CreateBitsrings()
    global M = []
    i = 1 
    while i <= n
        push!(M, CreateBitstring(i))
        i += 1
    end
end

#=

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

=#

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

function dom_vars()
    vars = []
    for p in dom_pairs
        append!(vars, p)
        #= if !(abs(p[1]) in vars)
            append!(vars, p[1])
        end
        if !(abs(p[2]) in vars)
            append!(vars, p[2])
        end =#
    end
    #println(length(vars))
    unique!(vars)
    return(length(vars))
end

function domposs(x,y)
    sup1 = true
    sup2 = true
    i = 1
    while i <= 32 && (sup1 || sup2)
        if M[x][i] - M[y][i] == 1
            sup2 = false
        elseif M[x][i] - M[y][i] == -1
            sup1 = false
        end
        i += 1
    end
    return(sup1 || sup2)
end
            
#global L
#=tic()
CreateHeaps()
println("Creating Heaps took: ")
toc()
tic()

test = []
for i in 1 : n
    append!(test, i)
end 
L = [[1,test]] =#
#global runtime = 0
#global printer = 1
#l = length(L)
function main()
    outputs = true
    global start = time()
    #tic()
    CreateBitsrings()
    #toc()
    global x = 1
    while x <= n && (time() - start) <= 10
        # println(x)
        minlength = 1000000
        minindex = 0
        for i in findnz(A[:, x])[1]
            if length(findnz(A[i, :])[1]) < minlength 
                minlength = length(findnz(A[i, :])[1])
                minindex = i
            end
        end
        if minindex != 0
            for y in findnz(A[minindex, :])[1]
                if y != x && domposs(x,y)
                    compare(y,x)
                    compareneg(y,x)
                end
            end
        end
        x += 1
    end
    if time() - start > 10
        #println("couldn't complete")
        outputs = false
    end
    #println("done")
    unique!(dom_pairs)
    #println("done")
    #println(length(dom_pairs))
    res = dom_vars()
    #println("done")
    #println(outputs,length(dom_pairs),res)
    return(outputs,length(dom_pairs),res)
end

cd("C:/Users/jacob/Desktop/Testing")

directory = readdir()

global directoryindex = 4
global const directorylength = length(directory)

while directoryindex <= directorylength
    global runtime = 0
    
    println(directory[directoryindex])
    model = read_from_file(directory[directoryindex])
    data = lp_matrix_data(model);
    
    global dom_pairs = []
    global A = data.A
    global c = data.c
    global n = length(c)
    global m = floor(Int, length(A) / n)
    global dom_pairs = []
    global zbitstring = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    global start1 = time()

    result = main()
    
    global runtime = time() - start1

    a = 1000
    open("0000resultsscip10.txt","a") do io
        println(io, "#", directory[directoryindex], "#",runtime , "#", result[1] , "#",  result[2], "#", result[3], "#", n, "#", m, "#", length(findnz(A)[3]), "#" )
    end

    println(directoryindex)

    global directoryindex += 1
end


    
#=
main()




model = read_from_file("C:/Users/Jacob/Desktop/Testing/mik-250-20-75-4.mps")
data = lp_matrix_data(model);









    lenbefore = length(dom_pairs)
    i1 = 1
    templ = length(L[iter][2])
    while i1 < templ #&& runtime <= 1000
        #global start = time()
        i2 = i1 + 1
        while i2 <= templ
            # println("comparing ", L[i][2][i1], " with ", L[i][2][i2])
            compare(L[iter][2][i1], L[iter][2][i2])
            compareneg(L[iter][2][i1], L[iter][2][i2])
            i2 += 1
        end
        # println("that was all for ", L[i][2][i1])
        i1 += 1
        #global runtime += time() - start
        #=if runtime >= printer*10
            println(runtime)
            global printer += 1
        end=#
    end
    if length(dom_pairs) >= lenbefore
        println(length(dom_pairs))
    end
    #println(i/l, "%")
    global iter += 1
end
println("Runtime was: ")
toc()
println(" there are this many pairs: ", length(dom_pairs))
#println(dom_pairs)
            
=#









