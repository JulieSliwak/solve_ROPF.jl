include(joinpath("D:\\repo", "ComplexOPF.jl","src", "PowSysMod_body.jl"))
include("SDP_decomposition_functions.jl")
include("solve_SDP.jl")
include("run_ampl_minlp.jl")
include("B&BandB_maxk_fixingsome1.jl")

struct ROPF_infos
    instance_name::String
    matpower_instance_path::String
    output_instance_path::String
    decomposition::String
    output_decomposition_path::String
    generation::String
    generation_files_path::String
end

struct BB_infos
    search_strategy::String
    branch_strategy::String
    seuil_u::Float64
    seuil_l::Float64
end


function construct_dat_file_ROPF(instance_name, matpower_instance_path, output_path)
    typeofinput = MatpowerRTEROPFSimpleInput
    OPFpbs = load_OPFproblems(typeofinput, matpower_instance_path, flag)
    # Bulding optimization problem
    pb_global = build_globalpb!(OPFpbs)
    pb_global_real = pb_cplx2real(pb_global)
    export_to_dat(pb_global_real, output_path, filename="$(instance_name).dat")
    return
end


function generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)
    sp = read_sparsity_pattern(matpower_instance_path)
    nb_nodes[instance] = size(sp,1)
    nb_edges[instance] = (nnz(sp) - size(sp,1))/2
    L, order, nb_added_edges = chordal_ext_cholesky(sp)
    H = construct_graph_from_matrix(L)
    cliques_dict = cliques_max(H,order)
    GC = weighted_graph(cliques_dict)
    A = Prim_algo(GC)
    isdir(joinpath(output_decomposition_path, "cliquetree_cholesky")) || mkpath(joinpath(output_decomposition_path, "cliquetree_cholesky"))
    f = open(joinpath(output_decomposition_path, "cliquetree_cholesky", "$(instance_name)_sdp_cliquetree.txt"), "w")
        for edge in LightGraphs.edges(A)
            clique1 = src(edge)
            clique2 = dst(edge)
            @printf(f, "%10s %20s\n", "B$clique1", "B$clique2")
        end
    close(f)
    isdir(joinpath(output_decomposition_path, "blocks_cholesky")) || mkpath(joinpath(output_decomposition_path, "blocks_cholesky"))
    f = open(joinpath(output_decomposition_path, "blocks_cholesky", "$(instance_name)_sdp_blocks.txt"), "w")
    for (clique, nodes_list) in cliques_dict
        for node in nodes_list
            @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(node)_Re")
            @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(node)_Im")
        end
    end
    close(f)
end


function solve1(ROPF)
    LB_plus, stat_plus = solve_SDP(ROPF, "plus")
    LB_minus, stat_minus = solve_SDP(ROPF, "minus")
    UB_minus = UB_plus = ""
    if !(stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT) && ! (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        println("$instance $generation : SDP relaxation probably infeasible ")
    end
    if (stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        UB_plus = solve_minlp(ROPF, "plus", [], Dict{String,Float64}())
    end
    if (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        UB_minus = solve_minlp(ROPF, "minus", [], Dict{String,Float64}())
    end
    println("Increase in generation : UB=$UB_plus ; LB = $LB_plus \n")
    println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
    return UB_plus, LB_plus, UB_minus, LB_minus
end


function solve2(ROPF, max_time)
    LB_plus, stat_plus = solve_SDP(ROPF, "plus")
    LB_minus, stat_minus = solve_SDP(ROPF, "minus")
    UB_minus = UB_plus = ""
    if !(stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT) && ! (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        println("$instance $generation : SDP relaxation probably infeasible ")
    end
    if (stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT)
        #B&B algo
        BB_parameters = BB_infos("deepfirst", "1", 0.9, 10^(-4))
        (UB_plus, nb_nodes, open_nodes) = BandB_maxk_fixingsome1and0(ROPF, "plus", BB_parameters, max_time)
    end
    if (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        #B&B algo
        BB_parameters = BB_infos("deepfirst", "1", 0.9, 10^(-4))
        (UB_minus, nb_nodes, open_nodes) = BandB_maxk_fixingsome1and0(ROPF, "minus", BB_parameters, max_time)
    end
    println("Increase in generation : UB=$UB_plus ; LB = $LB_plus \n")
    println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
    return UB_plus, LB_plus, UB_minus, LB_minus
end
