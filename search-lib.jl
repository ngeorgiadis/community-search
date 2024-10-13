include("domination-v2.jl")

using ProgressMeter
using CSV
using DataFrames
using Dates
using GraphIO.GML
using GraphPlot
using Graphs
using MD5
using MetaGraphs
using Statistics
using TOML
using JSON
using ArgParse


mutable struct HopSearchOptions
    hops::Int8
    HopSearchOptions() = new(1)
    HopSearchOptions(x) = new(x)
end

mutable struct RandomwalkSearchOptions
    steps::Int8
    iter::Int16
    RandomwalkSearchOptions() = new(15, 10)
    RandomwalkSearchOptions(x, y) = new(x, y)
end

mutable struct Config
    nodes_file::String
    edges_file::String
    query::Int64
    search_mode::String
    hop_search::HopSearchOptions
    randomwalk_search::RandomwalkSearchOptions

    # global pre computed domination scores file
    domination_scores_file::String

    Config(nodes_file, edges_file, query, search_mode) = new(
        nodes_file,
        edges_file,
        query,
        search_mode,
        HopSearchOptions(),
        RandomwalkSearchOptions(),
        ""
    )
end

function get_config()
    c = TOML.parsefile("config.toml")

    parsed = Config

    parsed = Config(
        c["nodes_file"],
        c["edges_file"],
        c["query"],
        "HOP_SEARCH",
    )

    if haskey(c, "search_mode")
        if c["search_mode"] == "HOP_SEARCH" || c["search_mode"] == "RANDOMWALK_SEARCH"
            parsed.search_mode = c["search_mode"]

            if c["search_mode"] == "HOP_SEARCH"
                parsed.hop_search.hops = c["hop_search"]["hops"]
            elseif c["search_mode"] == "RANDOMWALK_SEARCH"
                parsed.randomwalk_search.steps = c["randomwalk_search"]["steps"]
                parsed.randomwalk_search.iter = c["randomwalk_search"]["iter"]
            end
        else
            error("config value of search_mode can be 'HOP_SEARCH' OR 'RANDOMWALK_SEARCH' ")
        end
    end

    if haskey(c, "domination_scores_file")
        parsed.domination_scores_file = c["domination_scores_file"]
    end

    return parsed
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--query"
        help = "query ID"
        arg_type = Int
        default = -1
        "--hops"
        help = "number of hops in HOP_SEARCH mode"
        arg_type = Int
        default = -1
        "--mode"
        help = "search mode (HOP_SEARCH or RANDOMWALK_SEARCH)"
        arg_type = String
        default = ""
        "--steps"
        help = "steps of RANDOMWALK_SEARCH"
        arg_type = Int
        default = -1
        "--iter"
        help = "num of iterations in RANDOMWALK_SEARCH"
        arg_type = Int
        default = -1
    end
    return parse_args(s)
end

function find_communities(g, dsa, top, hops, check_points)
    results = []
    visited = Dict{Int64,Bool}()
    i = 1
    max_dom = dsa[1][:dom]

    check1 = Base.time()
    checpoint_times = []
    cpi = 1
    egotime = 0
    kcoretime = 0

    for (idx, n) in enumerate(dsa[1:top])

        if (idx >= check_points[cpi])
            p1 = Base.time() - check1
            push!(checpoint_times, p1)
            println("")
            println("checkpoint: $(sum(checpoint_times)), top-$(check_points[cpi]), egonet time: $(egotime), k-core time:$(kcoretime) ( $(length(results)) ), $(hops)")
            cpi += 1
            check1 = Base.time()
        end

        if haskey(visited, n[:id])
            continue
        end

        t0 = Base.time()

        # find egonet
        e1 = egonet(g, n[:id], hops)
        # end

        et = Base.time() - t0
        egotime = egotime + et

        v1 = filter(
            v -> !haskey(visited, get_prop(e1, v, :id)),
            vertices(e1),
        )

        e2 = e1[v1]
        if nv(e2) <= 1
            continue
        end


        t0 = Base.time()

        # find max k-core
        corenum = core_number(e2)
        k = maximum(corenum)
        max_k_core = findall(x -> x >= k, corenum)
        k1, _ = induced_subgraph(e2, max_k_core)
        # end

        kt = Base.time() - t0
        kcoretime = kcoretime + kt

        k2 = map(v -> props(k1, v), vertices(k1))

        for v in k2
            visited[v[:id]] = true
        end

        # get graph stats
        stats = get_graph_stats(k2, max_dom, k)
        stats["init"] = n[:id]
        stats["original_index"] = idx
        stats["egotime"] = et
        stats["kcoretime"] = kt

        if DEBUG
            stats["avg_clustering"] = mean(local_clustering_coefficient(k1))
            stats["density"] = density(k1)
            stats["avg_degree"] = 2 * ne(k1) / nv(k1)

            nodes = ""
            for n in k2
                nodes = nodes * "$(n[:id]), "
            end
            stats["nodes"] = nodes
        end

        push!(results, stats)

        if i % 1000 == 0
            print(".")
        end
        i += 1
    end

    return results, egotime, kcoretime
end

function create_graph(file)
    G = MetaGraph(1712433)

    open(file) do f
        for ln in eachline(f)
            p = split(ln, "\t")
            if length(p) == 3
                index1 = parse(Int64, SubString(p[1], 2))
                index2 = parse(Int64, p[2])
                add_edge!(G, index1, index2)
            end
        end
    end

    for (i, v) in enumerate(vertices(G))
        set_prop!(G, v, :id, v)
    end

    return G
end

function a_dominates_b(a, b)

    if size(a, 1) != size(b, 1)
        error("vectors should be the same length")
    end

    at_least_one_better = false
    for i in axes(a, 1)
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

"""
    domination_score(v, m)

    v is a Vector of length N
    m is a Matrix of size M x N

    Gets the domination score of a vector compared to other vectors
    passed as a matrix. This implementation is just a double loop 
    which is done for testing purposes and in case of large matrix will
    take long time to complete.
"""
function domination_score(v, m)
    score = 0
    for i in axes(m, 1)
        if a_dominates_b(v, m[i, 2:5])
            score = score + 1
        end
    end
    return score
end

function find_community(g, idx, hops)
    e = egonet(g, idx, hops)
    corenum = core_number(e)
    k = maximum(corenum)
    max_k_core = findall(x -> x >= k, corenum)
    k1, _ = induced_subgraph(e, max_k_core)
    return k1, k
end

function get_graph_stats(community, max_dom, max_core_number)

    stats = Dict{String,Any}()

    nodes_ds = map(x -> get_prop(community, x, :ds), vertices(community))

    N = length(nodes_ds)
    square_sum = 0

    for n in nodes_ds
        square_sum += (n - max_dom)^2
    end

    ratio_max_core = max_core_number / N
    stddev = sqrt(square_sum / N)

    stats["number_of_nodes"] = N
    stats["ratio_max_k_core"] = ratio_max_core
    stats["max_k_core"] = max_core_number
    stats["distance_from_max_ds"] = stddev
    stats["density"] = density(community)
    return stats
end

function graph_to_dict(g, rows, ds)
    nodes = map(x -> Dict{String,Any}(
            "id" => "node$(x)",
            "label" => rows[get_prop(g, x, :id)].name,
            "data" => Dict{String,Any}(
                "name" => rows[get_prop(g, x, :id)].name,
                "ds" => ds[get_prop(g, x, :id)],
                "init" => get_prop(g, x, :id),
            )), vertices(g))

    e1 = map(x -> Dict{String,Any}(
            "id" => "edge$(x[1])",
            "source" => "node$(src(x[2]))",
            "target" => "node$(dst(x[2]))",
        ), enumerate(edges(g)))

    res = Dict{String,Any}(
        "nodes" => nodes,
        "edges" => e1,
    )

    return res
end

function get_json_string(g, rows, ds)
    return JSON.json(graph_to_dict(g, rows, ds))
end

function get_graph_nodes_set(g)
    return [i for i in Set(map(x -> get_prop(g, x, :id), vertices(g)))]
end

function hop_search(g, id, hops, verbose=true)
    # find node index with id
    verbose && println("finding query index...")
    qIdx = findfirst(x -> get_prop(g, x, :id) == id, vertices(g))

    verbose && println("create egonet based on the query (node_id: $(id) index: $(qIdx)) with $(hops) hops...")
    return egonet(g, qIdx, hops)
end

function randomwalk_search(g, query, steps, iter)
    s = Vector{Int}()
    for i in 1:iter
        s = union(s, randomwalk(g, query, steps))
    end
    g1, _ = induced_subgraph(g, s)
    return g1
end

function find_community_randomwalk(g, idx, steps, iter)
    g1 = randomwalk_search(g, idx, steps, iter)

    corenum = core_number(g1)
    k = maximum(corenum)
    max_k_core = findall(x -> x >= k, corenum)
    k1, _ = induced_subgraph(g1, max_k_core)

    return k1, k
end