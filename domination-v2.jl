using TOML
using Dates
using DataFrames
using CSV
using Statistics

struct DatasetRow
    id::Int64
    name::String
    attrs::Vector{Int64}
end

mutable struct DatasetStats
    count::Int64
    max::Vector{Int64}
    min::Vector{Int64}
    histogram::Vector{Dict{Int64,Int64}}
end

struct DatasetPoint
    attrs::Vector{Int64}
    count::Int64
end

function read_dataset(csv_file::String)

    df = CSV.read(csv_file, DataFrame)
    df[!, :pi] = trunc.(Int64, df[!, :pi])

    return read_dataset(df)
end

function read_dataset(df::DataFrame)
    res = Dict{Int64,DatasetRow}()
    stats = DatasetStats(
        0,
        Int64[typemin(Int64), typemin(Int64), typemin(Int64), typemin(Int64)],
        Int64[typemax(Int64), typemax(Int64), typemax(Int64), typemax(Int64)],
        Dict{Int64,Int64}[
            Dict{Int64,Int64}(),
            Dict{Int64,Int64}(),
            Dict{Int64,Int64}(),
            Dict{Int64,Int64}(),
        ]
    )
    uniquePoints = Dict{String,Int64}()

    for (i, row) in enumerate(eachrow(df))
        id = row[:Id]
        name = row[:Label]
        attrs = Vector{Int64}(row[[:pc, :cn, :hi, :pi]])

        attrKey = ""
        for (i, a) in enumerate(attrs)

            if a > stats.max[i]
                stats.max[i] = a
            end

            if a < stats.min[i]
                stats.min[i] = a
            end

            if !haskey(stats.histogram[i], a)
                stats.histogram[i][a] = 0
            end
            stats.histogram[i][a] += 1

            attrKey *= "$(attrs[i])|"
        end

        attrKey = rstrip(attrKey, '|')
        if !haskey(uniquePoints, attrKey)
            uniquePoints[attrKey] = 0
        end
        uniquePoints[attrKey] += 1

        res[id] = DatasetRow(id, name, attrs)
    end

    stats.count = length(res)

    points = DatasetPoint[]

    for (k, v) in uniquePoints
        a = Int64[-1, -1, -1, -1]
        attrs = split(k, "|")
        for (i, v) in enumerate(attrs)
            a[i] = parse(Int64, v)
        end

        push!(points, DatasetPoint(a, v))
    end

    return res, stats, points
end

function datasetpointSortFn(x::DatasetPoint, y::DatasetPoint)
    s1 = sum(x.attrs)
    s2 = sum(y.attrs)
    if s1 > s2
        return true
    elseif s1 == s2
        sd1 = std(x.attrs)
        sd2 = std(y.attrs)
        return sd1 < sd2
    else
        return false
    end
end

function translateDatasetPoint(p::Vector{Int64}, stats::DatasetStats, gridSize::Vector{Int64})
    res = Vector{Int64}(undef, length(p))
    for (i, v) in enumerate(p)
        steps = (stats.max[i] - stats.min[i]) / gridSize[i]
        res[i] = 1 + trunc(Int64, v / steps)
    end

    return res
end

function getCoordinatekey(v::Vector{Int64})
    res = ""
    for i in eachindex(v)
        res *= "$(v[i])|"
    end
    return rstrip(res, '|')
end

function a_less_b(a::Vector{Int64}, b::Vector{Int64})
    for i in eachindex(a)
        if a[i] >= b[i]
            return false
        end
    end
    return true
end

function a_less_or_equal_b(a::Vector{Int64}, b::Vector{Int64})
    for i in eachindex(a)
        if a[i] > b[i]
            return false
        end
    end
    return true
end

function a_dominates_b(a, b)
    at_least_one_better = false
    for i in eachindex(a)
        if a[i] < b[i]
            return false
        end

        # it will come to this point 
        # if a[i] >= b[i]

        # mark if there is at least one element
        # a strictly better than b
        if a[i] > b[i]
            at_least_one_better = true
        end
    end
    return at_least_one_better
end

function calc_domination_score(stats::DatasetStats, points::Vector{DatasetPoint}, verbose=true)

    t0 = Base.time()

    sort!(points, lt=datasetpointSortFn)

    domination = Dict{String,Int64}()
    grid = Dict{String,Vector{DatasetPoint}}()

    gridSize = Int64[25, 25, 25, 25]

    for i in eachindex(points)
        v = points[i]
        coordinates = translateDatasetPoint(v.attrs, stats, gridSize)
        key = getCoordinatekey(coordinates)
        if !haskey(grid, key)
            grid[key] = DatasetPoint[]
        end
        push!(grid[key], v)
    end

    gridCoords = DatasetPoint[]
    for (k, v) in grid
        a = Int64[-1, -1, -1, -1]
        attrs = split(k, "|")
        for (i, v) in enumerate(attrs)
            a[i] = parse(Int64, v)
        end

        push!(gridCoords, DatasetPoint(a, 0))
    end

    sort!(gridCoords, lt=datasetpointSortFn)


    verbose && println("total cells to explore: $(length(gridCoords))")

    # main loop
    for i in eachindex(gridCoords)

        c = gridCoords[i]
        iKey = getCoordinatekey(c.attrs)
        point = c.attrs

        sumPoint = sum(point)

        baseScore = 0
        later = DatasetPoint[]

        for (_, j) in enumerate(gridCoords[i:end])
            if sum(j.attrs) > sumPoint
                continue
            end

            jk = getCoordinatekey(j.attrs)
            pointToCompareWith = j.attrs

            if a_less_b(pointToCompareWith, point)

                for (_, v) in enumerate(grid[jk])
                    baseScore += v.count
                end
            elseif a_less_or_equal_b(pointToCompareWith, point)
                push!(later, grid[jk]...)
            end
        end

        for (_, n) in enumerate(grid[iKey])

            nodeScore = baseScore

            for (_, l) in enumerate(later)
                if a_dominates_b(n.attrs, l.attrs)
                    nodeScore += l.count
                end
            end

            domination[getCoordinatekey(n.attrs)] = nodeScore
        end

        if i % 100 == 0
            verbose && print(".")
        end
    end

    verbose && println("\n\nTIME: domination score found in ", Base.time() - t0)

    return domination
end



