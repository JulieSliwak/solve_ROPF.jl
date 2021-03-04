function solve_minlp(ROPF, flag, fixing, index_var)
    instance = ROPF.instance_name
    generation = ROPF.generation
    root = pwd()
    # date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    if length(fixing) > 0
        outlog = joinpath(pwd(), "knitro_runs", "BB_$(instance)_$(generation)_$(flag).log")
    else
        outlog = joinpath(pwd(), "knitro_runs", "$(instance)_$(generation)_$(flag).log")
    end
    pb_path = ROPF.output_instance_path
    src_ampl_path = joinpath(pwd(), "src_ampl")

    cd(pb_path)

    instance_dat_file = "$(instance).dat"
    cp(instance_dat_file, "minlp_instance.dat", force=true)
    Pinput_lines = readlines(joinpath(ROPF.generation_files_path, "$(instance)_$(generation).csv"))
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

    f = open("fixing.dat", "w")
    all_fixing = true
    for (var, index) in index_var
        value = fixing[index]
        if value != -1
            write(f, "$var    $value \n")
        else
            all_fixing = false
        end
    end
    close(f)

    if all_fixing
        open("phase3.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "phase3.run"));")
        end
    else
        open("minlp.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "minlp.run"));")
        end
    end

    open("minlp.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp.mod"));")
    end

    try
        if all_fixing
            run(`cmd /c ampl phase3.run '>' $(outlog)`)
        else
            run(`cmd /c ampl minlp.run '>' $(outlog)`)
        end
        # run(`cmd /c ampl real_minlp.run `)
    catch
        @warn("AMPL/Knitro failed, returning.")
    end

    mv("knitro_solution.csv", joinpath(pwd(), "knitro_solution.csv"), force=true)
    # mv("knitro_solution.csv", "D:\\repo\\ROPF.jl\\knitro_optimal_solutions\\$(instance).csv", force=true)
    rm("minlp_instance.dat")
    rm("Pinput.dat")
    rm("lambda.dat")
    rm("fixing.dat")
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
