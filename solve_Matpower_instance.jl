include("generic_functions.jl")
#example
instance_name = "case1354pegase"
matpower_instance_path = "D:\\repo\\data\\data_Matpower\\matpower\\$(instance_name).m"
output_instance_path = "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"
output_decomposition_path = "D:\\repo\\data\\data_sdp"
generation = "generationrandom1"
generation_files_path = "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"
max_time = 120 #1 hour


#construct_dat_file_ROPF(instance_name, matpower_instance_path, output_instance_path)
#generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)
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
UB_plus, LB_plus, UB_minus, LB_minus = solve2(ROPF, max_time)
println("After B&B")
println("Increase in generation : UB=$UB_plus ; LB = $LB_plus")
println("Decrease in generation : UB=$UB_minus ; LB = $LB_minus \n")
