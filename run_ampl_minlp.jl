function run_knitro(pb_path::String, instance::String, src_ampl_path::String, flag::String, generation::String)
    root = pwd()
    # date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    outlog = "D:\\repo\\ROPF.jl\\knitro_runs\\Knitro_minlp_$(instance)_$(generation)_$(flag).log"

    cd(pb_path)

    instance_dat_file = "$(instance).dat"
    cp(instance_dat_file, "minlp_instance.dat", force=true)
    Pinput_lines = readlines("$(instance)_$(generation).csv")
    f = open("Pinput.dat", "w")
    for line in Pinput_lines
        temp = split(line, ';')
        var = temp[1]
        value = temp[2]
        Pmin = temp[4]
        Pmax = temp[5]
        write(f, "$var  $value  $Pmin   $Pmax \n")
    end
    close(f)
    f = open("lambda.dat", "w")
    write(f, "lambda_$flag   0.5")
    close(f)

    open("minlp.run", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp.run"));")
    end

    open("minlp.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp.mod"));")
    end

    try
        run(`cmd /c ampl minlp.run '>' $(outlog)`)
        # run(`cmd /c ampl real_minlp.run `)
    catch
        @warn("AMPL/Knitro failed, returning.")
    end

    #mv("knitro_solution.csv", "D:\\repo\\RTE_ROPFu\\knitro_solutions\\$(instance)_$(generation)_$(flag).csv", force=true)
    # mv("knitro_solution.csv", "D:\\repo\\RTE_ROPFu\\knitro_optimal_solutions\\$(instance).csv", force=true)
    rm("minlp_instance.dat")
    rm("Pinput.dat")
    rm("lambda.dat")
    cd(root)
    objective_value = 0.0
    status = ""
    time = 0.0
    lines = readlines(outlog)
    for line in lines
        splitted_line = split(line, ":")
        if splitted_line[1] == "EXIT"
            status = splitted_line[2]
        end
        splitted_line = split(line, "=")
        if splitted_line[1] == "Final objective value               "
            objective_value = parse(Float64,splitted_line[end])
        elseif splitted_line[1] == "Total program time (secs)           "
            t = split(splitted_line[2])[1]
            time += parse(Float64,t)
        end
    end

    if status == " Locally optimal solution found."
        return objective_value
    else
        return +Inf
    end
end
#
# pb_path = "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"
# # instances = ["case_ACTIVSg200", "case300","case_ACTIVSg500", "case1354pegase", "case1888rte", "case1951rte", "case_ACTIVSg2000",
# # "case2383wp", "case2736sp", "case2737sop", "case2746wop", "case2746wp", "case2848rte", "case2868rte",
# # "case2869pegase", "case3012wp", "case3120sp", "case3375wp" , "case6468rte","case6470rte",
# # "case6495rte", "case6515rte"#=,  "case9241pegase", "case13659pegase"=#]
# src_ampl_path = "D:\\repo\\RTE_ROPFu\\src_ampl"
# # generation = "generation1"
# # flag = "minus"
# instances = ["case2383wp"]
# generation = "generationrandom2"
# for instance in instances
#     obj_plus = run_knitro(pb_path, instance, src_ampl_path, flag, generation)
#     println(obj_plus)
# end
