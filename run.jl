using JSON3
include("search-lib.jl")

function read_test_authors()
    auth = JSON3.read("random_selected_nodes_100.json")
    auth_array = []
    for (k, v) in auth
        push!(auth_array, v)
    end
    ids = map(x -> x["id"], auth_array)
    return ids
end

function hop_search_wrapper(g, df, global_ds_indx, id, hops)

    g1 = hop_search(g, id, hops + 1, false)
    m = map(x -> get_prop(g1, x, :id), vertices(g1))

    # !!!
    # df needs to be unsorted
    # because the index is the id
    df_ego = df[m, :]
    rows, stats, points = read_dataset(df_ego)

    # println("calculating domination scores...")
    domination = calc_domination_score(stats, points, false)
    df_ego[!, :ds] .= -1

    ds = Dict{Int64,Int64}()
    for (k, v) in rows
        ds[k] = domination[getCoordinatekey(v.attrs)]
    end

    transform!(df_ego, [:Id] => ((x) -> map(y -> ds[y], x)) => :ds)
    # sort by domination score desc (higher at top)
    sort!(df_ego, [:ds], rev=true)
    max_ds = df_ego[1, :ds]

    # add the domination score in the graphs attribures
    for v in vertices(g1)
        set_prop!(g1, v, :ds, ds[get_prop(g1, v, :id)])
    end
    # println("searching for subgraphs in $(size(df_ego))")
    # starting from the highest domination score
    # create a subgraph (egonet) starting from ith node with distance c.hops

    id2idx = Dict{Int64,Int64}()
    for v in vertices(g1)
        id2idx[get_prop(g1, v, :id)] = v
    end
    res = Any[]
    i = 1
    lk = ReentrantLock()

    # println("evaluating each community in $(Threads.nthreads()) threads")
    Threads.@threads for row in eachrow(df_ego)
        rowid = row[:Id]

        if !haskey(id2idx, rowid)
            error("cannot find index for id: $(rowid)")
        end

        idx = id2idx[rowid]
        community, max_k_core = find_community(g1, idx, hops)

        qIdx = findfirst(x -> get_prop(community, x, :id) == id, vertices(community))
        if isnothing(qIdx)
            continue
        end

        g1_stats = get_graph_stats(community, max_ds, max_k_core)
        g1_stats["init"] = rowid
        g1_stats["idx"] = idx
        # g1_stats["json"] = get_json_string(community, rows, ds)
        g1_stats["nodes_set"] = get_graph_nodes_set(community)

        cds = map(x -> global_ds_indx[get_prop(community, x, :id)], vertices(community))
        g1_stats["avg_ds"] = sum(cds) / length(cds)

        lock(lk) do
            push!(res, g1_stats)
        end
    end
    return res
end

function random_walk_search_wrapper(g, df, global_ds_indx, id, len, iter)

    g1 = randomwalk_search(g, id, len, iter)
    m = map(x -> get_prop(g1, x, :id), vertices(g1))

    df_ego = df[m, :]
    rows, stats, points = read_dataset(df_ego)

    # println("calculating domination scores...")
    domination = calc_domination_score(stats, points, false)
    df_ego[!, :ds] .= -1

    ds = Dict{Int64,Int64}()
    for (k, v) in rows
        ds[k] = domination[getCoordinatekey(v.attrs)]
    end

    transform!(df_ego, [:Id] => ((x) -> map(y -> ds[y], x)) => :ds)
    # sort by domination score desc (higher at top)
    sort!(df_ego, [:ds], rev=true)
    max_ds = df_ego[1, :ds]

    # add the domination score in the graphs attribures
    for v in vertices(g1)
        set_prop!(g1, v, :ds, ds[get_prop(g1, v, :id)])
    end
    # println("searching for subgraphs in $(size(df_ego))")
    # starting from the highest domination score
    # create a subgraph (egonet) starting from ith node with distance c.hops

    id2idx = Dict{Int64,Int64}()
    for v in vertices(g1)
        id2idx[get_prop(g1, v, :id)] = v
    end
    res = Any[]
    i = 1
    lk = ReentrantLock()

    # println("evaluating each community in $(Threads.nthreads()) threads")
    Threads.@threads for row in eachrow(df_ego)
        #
        rowid = row[:Id]
        if !haskey(id2idx, rowid)
            error("cannot find index for id: $(rowid)")
        end

        idx = id2idx[rowid]
        community, max_k_core = find_community_randomwalk(
            g1,
            idx,
            len,
            iter
        )

        qIdx = findfirst(x -> get_prop(community, x, :id) == id, vertices(community))
        if isnothing(qIdx)
            continue
        end

        g1_stats = get_graph_stats(community, max_ds, max_k_core)
        g1_stats["init"] = rowid
        g1_stats["idx"] = idx
        # g1_stats["json"] = get_json_string(community, rows, ds)
        g1_stats["nodes_set"] = get_graph_nodes_set(community)

        cds = map(x -> global_ds_indx[get_prop(community, x, :id)], vertices(community))
        g1_stats["avg_ds"] = sum(cds) / length(cds)

        lock(lk) do
            push!(res, g1_stats)
        end
    end
    return res
end

function store_results(res, name)
    if length(res) > 0
        try
            df_res = DataFrame(res)
            # println("saving csv file... dataframe size: $(size(df_res))")

            sort!(df_res, [:distance_from_max_ds])

            open("out/run-$(name).csv", "w") do output
                CSV.write(output, df_res, delim=";")
            end
        catch
            println("saving: $(name) failed")
        end
    else
        println("no results [$(name)]")

    end
end

function main()
    ids = read_test_authors()
    println("starting script...")
    println("$(Threads.nthreads()) threads available")


    println("reading config...")
    c = get_config()

    println("loading graph...")
    g = @time create_graph(c.edges_file)

    print("loading nodes...")
    df = CSV.read(c.nodes_file, DataFrame)
    df[!, :pi] = trunc.(Int64, df[!, :pi])
    println("done!")

    println("loading global domination scores...")
    global_ds_indx = Dict{Int,Int}()
    if c.domination_scores_file != ""
        ds_df = CSV.read(c.domination_scores_file, DataFrame, delim="\t", header=false)
        global_ds_indx = Dict{Int,Int}([r[1] => r[2] for r in eachrow(ds_df)])
    end


    # hop search

    times_one_hop = DataFrame([[], [], []], ["id", "runtime", "communities"])
    times_two_hops = DataFrame([[], [], []], ["id", "runtime", "communities"])

    # 1 hop search
    println("1 hop")
    for id in ids
        # run 
        # measure
        elapsedTime = @elapsed res = hop_search_wrapper(g, df, global_ds_indx, parse(Int64, id), 1)
        println("$(id) -> $(elapsedTime)")

        # store results
        fname = "$(id)-hop-1"
        store_results(res, fname)
        push!(times_one_hop, [id, elapsedTime, length(res)])
    end

    # 2 hop search
    println("2 hop")
    for id in ids
        # run 
        # measure
        elapsedTime = @elapsed res = hop_search_wrapper(g, df, global_ds_indx, parse(Int64, id), 2)
        println("$(id) -> $(elapsedTime)")

        # store results
        fname = "$(id)-hop-2"

        store_results(res, fname)
        push!(times_two_hops, [id, elapsedTime, length(res)])
    end


    open("out/run-times-one-hop.csv", "w") do output
        CSV.write(output, times_one_hop, delim=";")
    end

    open("out/run-times-two-hops.csv", "w") do output
        CSV.write(output, times_two_hops, delim=";")
    end

    #random walk search

    paths = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    iters = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    for p in paths
        for i in iters

            times = DataFrame([[], [], []], ["id", "runtime", "communities"])

            path_len = p
            iter = i
            println("RW $(path_len) $(iter)")
            for id in ids
                # run 
                # measure
                elapsedTime = @elapsed res = random_walk_search_wrapper(g, df, global_ds_indx, parse(Int64, id), path_len, iter)
                println("$(id) -> $(elapsedTime)")

                # store results
                fname = "$(id)-walk-$(path_len)-$(iter)"
                store_results(res, fname)
                push!(times, [id, elapsedTime, length(res)])
            end

            open("out/run-times-walk-$(path_len)-$(iter).csv", "w") do output
                CSV.write(output, times, delim=";")
            end
        end
    end
end

@time main()
