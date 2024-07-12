using JuMP
using HiGHS
using SparseArrays
using Plots

f = open("C:/Users/Jacob/Desktop/Testing/0000resultsscip10.txt", "r")
line_count = 0
Data = []
for lines in readlines(f)
    global line_count += 1
    push!(Data, [lines])
    #println(lines)
    global line_count += 1
end

close(f)

f = open("C:/Users/Jacob/Desktop/Testing/0000resultsmine10.txt", "r")
line_count = 0
Data2 = []
for lines in readlines(f)
    global line_count += 1
    push!(Data2, [lines])
    #println(lines)
    global line_count += 1
end

close(f)

global BetterData = []
global geometricmean1 = 1.0
global geometricmean2 = 1.0
global i = 1
while i <= 240
    println(i)
    temp = []
    pos = 2
    while Data[i][1][pos] != '#'
        #println(Data[i][pos])
        pos += 1
    end
    push!(temp, Data[i][1][2:pos-1])
    pos2 = pos + 1
    while Data[i][1][pos2] != '#'
        pos2 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos+1:pos2-1])])
    pos2 += 1
    pos3 = pos2 +1
    while Data[i][1][pos3] != '#'
        pos3 += 1
    end
    push!(temp, [parse(Bool, Data[i][1][pos2:pos3-1])])
    pos3 += 1
    pos4 = pos3 + 1
    while Data[i][1][pos4] != '#'
        pos4 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos3:pos4-1])])
    pos4 += 1
    pos5 =  pos4 + 1
    while Data[i][1][pos5] != '#'
        pos5 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos4:pos5-1])])
    pos5 += 1
    pos6 =  pos5 + 1
    while Data[i][1][pos6] != '#'
        pos6 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos5:pos6-1])])
    pos6 += 1
    pos7 =  pos6 + 1
    while Data[i][1][pos7] != '#'
        pos7 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos6:pos7-1])])
    pos7 += 1
    pos8 =  pos7 + 1
    while Data[i][1][pos8] != '#'
        pos8 += 1
    end
    push!(temp, [parse(Float64, Data[i][1][pos7:pos8-1])])
    global i += 1
    push!(BetterData, temp)
end

global i = 1
while i <= 240
    println(i)
    temp = []
    pos = 2
    while Data2[i][1][pos] != '#'
        #println(Data[i][pos])
        pos += 1
    end
    #push!(temp, Data[i][1][2:pos-1])
    pos2 = pos + 1
    while Data2[i][1][pos2] != '#'
        pos2 += 1
    end
    push!(BetterData[i], [parse(Float64, Data2[i][1][pos+1:pos2-1])])
    pos2 += 1
    pos3 = pos2 +1
    while Data2[i][1][pos3] != '#'
        pos3 += 1
    end
    push!(BetterData[i], [parse(Float64, Data2[i][1][pos2:pos3-1])])
    pos3 += 1
    pos4 = pos3 + 1
    while Data2[i][1][pos4] != '#'
        pos4 += 1
    end
    push!(BetterData[i], [parse(Float64, Data2[i][1][pos3:pos4-1])])
    pos4 += 1
    pos5 =  pos4 + 1
    while Data2[i][1][pos5] != '#'
        pos5 += 1
    end
    push!(BetterData[i], [parse(Float64, Data2[i][1][pos4:pos5-1])])
    global i +=1
end

function geometricmean()
    osol1 = 0
    osol2 = 0
    nsol = 0
    faster1 = 0
    i = 1
    j = 1
    while i <= length(BetterData)
        if BetterData[i][2][1] >= 9
            global geometricmean1 = geometricmean1^((j -1)/j)*((BetterData[i][2][1]+1))^(1/j)
            global geometricmean2 = geometricmean2^((j -1)/j)*((min(BetterData[i][9][1],10)+1))^(1/j)
            j += 1
        end
        i += 1
    end
    geometricmean3 = 1
    geometricmean4 = 1
    geometricmean5 = 1
    geometricmean6 = 1
    mybetter = 0
    scipbetter = 0
    i = 1
    j=1
    while i <= length(BetterData)
        if  BetterData[i][2][1] >= 9
            geometricmean3 = geometricmean3^((j -1)/j)*((BetterData[i][4][1]+1))^(1/j)
            geometricmean4 = geometricmean4^((j -1)/j)*((BetterData[i][5][1]+1))^(1/j)
            geometricmean5 = geometricmean5^((j -1)/j)*((BetterData[i][10][1]+1))^(1/j)
            geometricmean6 = geometricmean6^((j -1)/j)*((BetterData[i][11][1]+1))^(1/j)
            if BetterData[i][11][1] > BetterData[i][5][1]
                mybetter += 1
            elseif BetterData[i][11][1] < BetterData[i][5][1]
                scipbetter += 1
            end
            j += 1
        end
        i += 1
    end
    println("Ratio between Geometric Means is: ", geometricmean1/geometricmean2, " solved by both: ", j-1)
    println("gm of scip pairs: ", geometricmean3)
    println("gm of scip vars: ", geometricmean4)
    println("gm of my pairs: ", geometricmean5)
    println("gm of my vars: ", geometricmean6)
    println("my better vers: ", mybetter)
    println("scip better vars: ", scipbetter)
    i = 1
    totallenmy = 0
    totallenba = 0
    totalvarmy = 0
    totalvarba = 0
    tpm = 0
    tpb = 0
    while i <= length(BetterData)
        totallenmy += BetterData[i][9][1]
        totallenba += BetterData[i][2][1]
        totalvarmy += BetterData[i][11][1]
        totalvarba += BetterData[i][5][1]
        tpm += BetterData[i][10][1]
        tpb += BetterData[i][4][1]
        i += 1
    end
    println("my ratio of var per s ", totalvarmy/totallenmy)
    println("my ratio of pairs per s ", tpm/totallenmy)
    println("base ratio of var per s ", totalvarba/totallenba)
    println("base ratio of pairs per s ", tpb/totallenba)
    println(totallenba)
    println(totallenmy)

    #=
    sort!(BetterData, by = x -> x[7][1])
    i = floor(Int, 3*length(BetterData)/4)
    j = 1
    while i <= length(BetterData)
        if BetterData[i][3][1] < 1000 && BetterData[i][4][1] < 1000
            global geometricmean1 = geometricmean1^((j -1)/j)*((BetterData[i][3][1]+1))^(1/j)
            global geometricmean2 = geometricmean2^((j -1)/j)*((BetterData[i][4][1]+1))^(1/j)
            j += 1
        elseif  BetterData[i][3][1] >= 1000 && BetterData[i][4][1] >= 1000
            nsol +=1
        elseif BetterData[i][3][1] >= 1000
            osol2 += 1
        else 
            osol1 +=1
        end
        i += 1
    end
    println("Ratio between Geometric Means on larger instances is: ", geometricmean1/geometricmean2, " number of instances: ", j-1)
    sort!(BetterData, by = x -> x[5][1]*x[6][1]/x[7][1])
    i = floor(Int, 3*length(BetterData)/4)
    j = 1
    while i <= length(BetterData)
        if BetterData[i][3][1] < 1000 && BetterData[i][4][1] < 1000
            global geometricmean1 = geometricmean1^((j -1)/j)*((BetterData[i][3][1]+1))^(1/j)
            global geometricmean2 = geometricmean2^((j -1)/j)*((BetterData[i][4][1]+1))^(1/j)
            j += 1
        elseif  BetterData[i][3][1] >= 1000 && BetterData[i][4][1] >= 1000
            nsol +=1
        elseif BetterData[i][3][1] >= 1000
            osol2 += 1
        else 
            osol1 +=1
        end
        i += 1
    end
    println("Ratio between Geometric Means on denser instances is: ", geometricmean1/geometricmean2, " number of instances: ", j-1) =#
end

geometricmean()
global minrunningtime = 1000
for i in BetterData
    if i[3][1] > 0 && i[3][1] < minrunningtime
        global minrunningtime = i[3][1]
    end
    if i[4][1] > 0 && i[4][1] < minrunningtime
        global minrunningtime = i[3][1]
    end
end
global maxratio = 1
global minratio = 1
for i in BetterData
    if max(minrunningtime, i[3][1])/ max(minrunningtime, i[4][1]) > maxratio
        global maxratio = max(minrunningtime, i[3][1])/ max(minrunningtime, i[4][1])
    end
    if max(minrunningtime, i[4][1])/ max(minrunningtime, i[3][1]) > minratio
        global minratio = max(minrunningtime, i[4][1])/ max(minrunningtime, i[3][1])
    end
end
println("Maxratio is: ", maxratio)
println("minratio is: ", minratio)

function plotbynonzeros()
    A = []
    B = []
    i = 1
    while i <= length(BetterData)
        push!(A, BetterData[i][2][3])
        push!(B, BetterData[i][2][4])
        i += 1
    end
    scatter(B,A)
    scatter!(xscale=:log10)
end

function plotbydensity()
    A = []
    B = []
    i = 1
    while i <= length(BetterData)
        push!(A, BetterData[i][2][3])
        push!(B, (BetterData[i][2][5])^(-1))
        i += 1
    end
    scatter!(B,A)
end

function plotbyn()
    A = []
    B = []
    i = 1
    while i <= length(BetterData)
        push!(A, BetterData[i][2][3])
        push!(B, BetterData[i][2][6])
        i += 1
    end
    scatter(B,A)
    scatter!(xscale=:log10)
end

function plotbyasymptotic()
    A = []
    B = []
    i = 1
    while i <= length(BetterData)
        push!(A, BetterData[i][2][3])
        push!(B, 1 + min(BetterData[i][2][1],BetterData[i][2][2]))
        i += 1
    end
    scatter(B,A)
    scatter!(xscale=:log10)
end

function plotbothbynz()
    A = []
    B = []
    C = []
    i = 1
    while i <= length(BetterData)
        push!(A, BetterData[i][2][4])
        push!(B, BetterData[i][2][1]+ 1)
        push!(C, BetterData[i][2][2]+ 1)
        i += 1
    end
    scatter(A,[B C])
    scatter!(xscale=:log10, yscale=:log10)
end

#sort!(BetterData, by = x -> x[2][3])
#println(BetterData)
#println((geometricmean1 - 1)/(geometricmean2 - 1))
#println((geometricmean1 - 1))
#println((geometricmean2 - 1))
#plotbynonzeros()
#plotbydensity()
#plotbyn()
#plotbyasymptotic()