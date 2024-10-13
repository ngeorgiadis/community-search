include("search-lib.jl")



function search()
    println("starting script...")

    c = get_config()
    # println("config: ", c)
    mkpath("out")

    println("check if command line options exists...")
    parsed_args = parse_commandline()

    if parsed_args["mode"] != ""
        c.search_mode = parsed_args["mode"]
    end

    if parsed_args["query"] > 0
        c.query = parsed_args["query"]
    end

    if c.search_mode == "HOP_SEARCH"
        if parsed_args["hops"] > 0
            c.hop_search.hops = parsed_args["hops"]
        end
    elseif c.search_mode == "RANDOMWALK_SEARCH"
        if parsed_args["steps"] > 0
            c.randomwalk_search.steps = parsed_args["steps"]
        end
        if parsed_args["iter"] > 0
            c.randomwalk_search.iter = parsed_args["iter"]
        end
    end

    println("loading graph...")
    g = @time create_graph(c.edges_file)
    # println("done!")

    print("loading nodes...")
    df = CSV.read(c.nodes_file, DataFrame)
    df[!, :pi] = trunc.(Int64, df[!, :pi])
    println("done!")

    # println("create egonet based on the query (node_id: $(c.query)) with $(c.hops + 1) hops...")
    # g1 = @time egonet(g, c.query, c.hops + 1)

    g1 = Nothing
    if c.search_mode == "HOP_SEARCH"
        g1 = hop_search(g, c.query, c.hop_search.hops + 1)
    elseif c.search_mode == "RANDOMWALK_SEARCH"
        g1 = randomwalk_search(g, c.query, c.randomwalk_search.steps, c.randomwalk_search.iter)
    else
        error("no search mode selected")
    end

    println("g1 subgraph contains $(nv(g1)) vertices and $(ne(g1)) edges.")

    print("saving graph to file...")
    if c.search_mode == "HOP_SEARCH"
        savegraph("out/$(c.query)-egograph.dot", g1, MetaGraphs.DOTFormat())
    elseif c.search_mode == "RANDOMWALK_SEARCH"
        savegraph("out/$(c.query)-walk.dot", g1, MetaGraphs.DOTFormat())
    end
    println("done!")

    m = map(x -> get_prop(g1, x, :id), vertices(g1))
    df_ego = df[m, :]

    rows, stats, points = read_dataset(df_ego)

    println("calculating domination scores...")
    domination = @time calc_domination_score(stats, points)

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


    println("searching for subgraphs in $(size(df_ego))")
    # starting from the highest domination score
    # create a subgraph (egonet) starting from ith node with distance c.hops

    id2idx = Dict{Int64,Int64}()
    for v in vertices(g1)
        id2idx[get_prop(g1, v, :id)] = v
    end

    res = Any[]
    i = 1
    lk = ReentrantLock()

    println("evaluating each community in $(Threads.nthreads()) threads")
    for row in eachrow(df_ego)
        #
        id = row[:Id]
        if !haskey(id2idx, id)
            error("cannot find index for id: $(id)")
        end

        idx = id2idx[id]

        if c.search_mode == "HOP_SEARCH"
            community, max_k_core = find_community(g1, idx, c.hop_search.hops)
        else
            community, max_k_core = find_community_randomwalk(
                g1,
                idx,
                c.randomwalk_search.steps,
                c.randomwalk_search.iter)
        end

        qIdx = findfirst(x -> get_prop(community, x, :id) == c.query, vertices(community))
        if isnothing(qIdx)
            # println("subgraph starting from $(id) doesnt contain query node. skipping...")
            continue
        end

        g1_stats = get_graph_stats(community, max_ds, max_k_core)
        g1_stats["init"] = id
        g1_stats["idx"] = idx
        g1_stats["json"] = get_json_string(community, rows, ds)
        push!(res, g1_stats)

        # lock(lk) do
        #     if i % 100 == 0
        #         print(".")
        #     end
        #     i += 1
        # end
    end

    if length(res) > 0
        # stmp = Dates.format(now(), "yyyymmdd-HHMMSS")
        try
            df_res = DataFrame(res)
            println("saving csv file... dataframe size: $(size(df_res))")
            sort!(df_res, [:distance_from_max_ds])

            s = c.search_mode === "HOP_SEARCH" ? "" : "walk"
            open("out/q-$(c.query)-$( s === "" ? c.hop_search.hops : s).csv", "w") do output
                CSV.write(output, df_res, delim=";")
            end
        catch
            println("PROCESS FAILED")
        end
    else
        println("no results. please try with another set of parameters")
    end
end

function main()
    search()
end



@time "Total process time" main()