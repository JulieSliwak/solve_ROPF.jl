include(joinpath("D:\\repo", "ComplexOPF.jl","src", "PowSysMod_body.jl"))
include("SDP_decomposition_functions.jl")
include("solve_SDP.jl")
include("run_ampl_minlp.jl")
include("B&BandB_maxk_fixingsome1.jl")

struct ROPF_instance
    instance_name::String
    matpower_instance_path::String
    output_instance_path::String
    formulation::String
    output_decomposition_path::String
    generation::String
    generation_files_path::String
end


function construct_ROPF_instance(instance_name, matpower_instance_path, output_path)
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


function solve1(instance, output_instance_path, output_decomposition_path, generation_files_path, generation)
    formulation = "cholesky"
    flag = "plus"
    LB_plus, stat_plus = (output_instance_path, output_decomposition_path, generation_files_path, formulation, instance, generation, flag)
    flag = "minus"
    LB_minus, stat_minus = (output_instance_path, output_decomposition_path, generation_files_path, formulation, instance, generation, flag)
    UB_minus = UB_plus = ""
    if !(stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT) && ! (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        println("$instance $generation : SDP relaxation probably infeasible ")
    end
    if (stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        src_ampl_path = joinpath(pwd(), "src_ampl")
        flag = "plus"
        UB_plus = run_knitro(output_instance_path, instance, src_ampl_path, flag, generation)
    end
    if (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        src_ampl_path = joinpath(pwd(), "src_ampl")
        flag = "minus"
        UB_minus = run_knitro(output_instance_path, instance, src_ampl_path, flag, generation)
    end
    println("Increase in generation : UB=$UB_plus ; LB = $LB_plus \n")
    println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
    return UB_plus, LB_plus, UB_minus, LB_minus
end


function solve2(instance, output_instance_path, output_decomposition_path, generation, max_time)
    formulation = "cholesky"
    flag = "plus"
    LB_plus, stat_plus = (output_instance_path, output_decomposition_path, generation_files_path, formulation, instance, generation, flag)
    flag = "minus"
    LB_minus, stat_minus = (output_instance_path, output_decomposition_path, generation_files_path, formulation, instance, generation, flag)
    UB_minus = UB_plus = ""
    if !(stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT) && ! (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        println("$instance $generation : SDP relaxation probably infeasible ")
    end
    if (stat_plus == MOI.FEASIBLE_POINT || stat_plus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        src_ampl_path = joinpath(pwd(), "src_ampl")
        flag = "plus"
        search_strategy = "deepfirst"
        branch_strategy = "1"
        max_var_1 = Inf
        seuil = 0.9
        t = @elapsed (UB_plus, nb_nodes, open_nodes) = BandB_maxk_fixingsome1(output_decomposition_path, formulation, instance, generation, flag, search_strategy, branch_strategy, max_var_1, max_time, seuil)
        println("Time : $t")
    end
    if (stat_minus == MOI.FEASIBLE_POINT || stat_minus == MOI.NEARLY_FEASIBLE_POINT)
        #solve MINLP with Knitro
        src_ampl_path = joinpath(pwd(), "src_ampl")
        flag = "minus"
        search_strategy = "deepfirst"
        branch_strategy = "1"
        max_var_1 = Inf
        seuil = 0.9
        t = @elapsed (UB_minus, nb_nodes, open_nodes) = BandB_maxk_fixingsome1(output_decomposition_path, formulation, instance, generation, flag, search_strategy, branch_strategy, max_var_1, max_time, seuil)
        println("Time : $t")
    end
    println("Increase in generation : UB=$UB_plus ; LB = $LB_plus \n")
    println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
    return UB_plus, LB_plus, UB_minus, LB_minus
end
