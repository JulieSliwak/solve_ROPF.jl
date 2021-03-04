include("generic_functions.jl")
#example
instance_name = "case14"
matpower_instance_path = "D:\\repo\\data\\data_Matpower\\matpower"
output_instance_path = "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"
output_decomposition_path = "D:\\repo\\data\\data_sdp"
generation = "generation1"
generation_files_path = "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"
max_time = 3600 #1 hour



construct_ROPF_instance(instance_name, matpower_instance_path, output_instance_path)
generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)
ROPF = ROPF_infos(instance_name,
matpower_instance_path,
output_instance_path,
"cholesky",
output_decomposition_path,
generation,
generation_files_path)
UB_plus, LB_plus, UB_minus, LB_minus = solve1(ROPF)
UB_plus, LB_plus, UB_minus, LB_minus = solve2(ROPF, max_time)
