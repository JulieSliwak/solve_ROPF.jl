include("generic_functions.jl")
#example
instances = ["case14"]
output_instance_path = "..\\data_ROPF"
output_decomposition_path = "..\\data_sdp"
list_generation = ["generation1"]
generation_files_path = "..\\data_ROPF"
max_time = 600 #10 minutes

for instance_name in instances
    matpower_instance_path = "..\\matpower_data\\$(instance_name).m"
    for generation in list_generation
        construct_dat_file_ROPF(instance_name, matpower_instance_path, output_instance_path)
        generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)
        ROPF = ROPF_infos(instance_name,
        matpower_instance_path,
        output_instance_path,
        "cholesky",
        output_decomposition_path,
        generation,
        generation_files_path)
        UB_plus, LB_plus, UB_minus, LB_minus = solve1(ROPF)
        println("Increase in generation : UB=$UB_plus ; LB = $LB_plus")
        println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
        if UB_minus <= UB_plus
            cp("solutions\\solve1_solution_$(instance_name)_minus.csv", "solutions\\solve1_BEST_solution_$(instance_name).csv")
        else
            cp("solutions\\solve1_solution_$(instance_name)_plus.csv", "solutions\\solve1_BEST_solution_$(instance_name).csv")
        end
        UB_plus, LB_plus, UB_minus, LB_minus = solve2(ROPF, max_time)
        println("After B&B")
        println("Increase in generation : UB=$UB_plus ; LB = $LB_plus")
        println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
        if UB_minus <= UB_plus
            cp("solutions\\BB_solution_$(instance_name)_minus.csv", "solutions\\solve2_BEST_solution_$(instance_name).csv")
        else
            cp("solutions\\BB_solution_$(instance_name)_plus.csv", "solutions\\solve2_BEST_solution_$(instance_name).csv")
        end
    end
end
