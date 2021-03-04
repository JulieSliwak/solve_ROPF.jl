using JuMP, Mosek, MosekTools, LinearAlgebra, Xpress, Statistics

function read_dat_file(dat_file_path)
    data = readdlm(dat_file_path, comments=true)
    SDP_var_list = []
    Sgen_var_list = []
    Bin_var_list = []
    λ = ""
    i = 1
    while data[i,1] == "VAR_TYPE"
        varname = data[i,3]
        if data[i,2] == "REAL"
            if varname[1:6] == "lambda"
                λ = varname
            elseif varname[1:3] == "Gen"
                push!(Sgen_var_list, varname)
            else
                push!(SDP_var_list, varname)
            end
        elseif data[i,2] == "BOOL"
            push!(Bin_var_list, varname)
        end
        i+=1
    end
    dict_quad_ctr = Dict{String, Dict{Tuple{String,String}, Float64}}()
    dict_linear_ctr = Dict{String, Dict{String,Float64}}()
    dict_bounds_ctr = Dict{String, Dict{String,Float64}}()
    dict_constants_ctr = Dict{String, Float64}()
    dict_Bin_ctr = Dict{String, Dict{String, Float64}}()
    dict_MONO = Dict{String,Array{Any}}()
    for k in i:size(data,1)
        type_elem = data[k,1]
        elem_name = data[k,2]
        var1 = data[k,3]
        var2 = data[k,4]
        value = data[k,5]
        if type_elem == "QUAD"
            if !haskey(dict_quad_ctr, elem_name)
                dict_quad_ctr[elem_name] = Dict((var1,var2)=>value)
            else
                dict_quad_ctr[elem_name][(var1,var2)] = value
            end
        elseif type_elem == "LIN"
            if !haskey(dict_linear_ctr, elem_name)
                dict_linear_ctr[elem_name] = Dict(var2 => value)
            else
                dict_linear_ctr[elem_name][var2] = value
            end
        elseif type_elem == "MONO"
            if !haskey(dict_Bin_ctr, elem_name)
                dict_Bin_ctr[elem_name] = Dict(var1 => value)
            else
                dict_Bin_ctr[elem_name][var1] = value
            end
        elseif type_elem == "LB" || type_elem == "UB"
            if !haskey(dict_bounds_ctr, elem_name)
                dict_bounds_ctr[elem_name] = Dict(type_elem => value)
            else
                dict_bounds_ctr[elem_name][type_elem] = value
            end
        elseif type_elem == "CONST"
            dict_constants_ctr[elem_name] = value
        elseif type_elem == "MONO_DEF"
            if !haskey(dict_MONO, elem_name)
                if value == 1
                    dict_MONO[elem_name] = [var1]
                elseif value == 2
                    dict_MONO[elem_name] = [(var1,var1)]
                end
            else
                if value == 1
                    push!(dict_MONO[elem_name], var1)
                elseif value == 2
                    push!(dict_MONO[elem_name], (var1,var1))
                end
            end
        end
    end
    # println(dict_bounds_ctr)

    return λ, Sgen_var_list, SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, dict_linear_ctr

end

function construct_SDP(blocks_dict, CLIQUE_TREE, Pinput_csv_file, flag, Sgen_var_list, SDP_var_list, Bin_var_list,
   dict_quad_ctr, dict_linear_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, solution_file)

   lines = readlines(Pinput_csv_file)
   NB_BLOCKS = length(blocks_dict)
    #initialize
    # m = Model(with_optimizer(Mosek.Optimizer, LOG=0))
    m = Model(with_optimizer(Mosek.Optimizer))
    #variables
    λ_and_Sgen_jumpvar = Dict(varname => @variable(m, base_name="$varname") for varname in Sgen_var_list)
    if flag == "plus"
      λ = "lambda_plus"
      λ_and_Sgen_jumpvar[λ] = @variable(m, base_name="$λ", lower_bound = 0.5, upper_bound=1)
    elseif flag == "minus"
      λ = "lambda_minus"
      λ_and_Sgen_jumpvar[λ] = @variable(m, base_name="$λ", lower_bound = 0, upper_bound=0.5)
   else
     exit("Problem with flag")
   end
   constraints_ref=Dict{String,JuMP.ConstraintRef}()
   for line in lines
       temp = split(line, ';')
       gen_name = split(temp[1])[1]
       Pinput = parse(Float64,temp[2])
       empty = temp[3]
       Pmin = parse(Float64,temp[4])
       Pmax = parse(Float64,temp[5])
       if Pinput > Pmax
           Pinput = Pmax
       elseif Pinput < Pmin
           Pinput = Pmin
       end
       Pgen_var = λ_and_Sgen_jumpvar[gen_name]
       bus = split(gen_name, '_')[end-1]
       ctrname = "_$(bus)_Gen_Activepower_Re"
       if flag == "minus"
         λminus = λ_and_Sgen_jumpvar[λ]
         constraints_ref[ctrname] = @constraint(m, -Pgen_var+ (Pmin+2*(Pinput-Pmin)*λminus) == 0)
      elseif flag == "plus"
        λplus = λ_and_Sgen_jumpvar[λ]
        constraints_ref[ctrname] = @constraint(m, -Pgen_var + (2*Pinput-Pmax+2*(Pmax- Pinput)*λplus) ==0)
      end
   end


    jumpBinvar = Dict{String, VariableRef}()
    uvar = Dict{String, Tuple{VariableRef, JuMP.GenericAffExpr}}()
    for bin_var in Bin_var_list
      jumpBinvar[bin_var] = @variable(m, base_name = "ξ_$bin_var", lower_bound = 0)
    end

    jumpX = Dict{String,Array{VariableRef,2}}()
    coeff_block = Dict{Tuple{String,String}, Set{String}}()
    mat_var = Dict{Tuple{String,String}, Dict{String,Any}}()

    for (block, var_list) in blocks_dict
      for var1 in var_list
        for var2 in var_list
          mat_var[(var1,var2)] = Dict{String,Any}()
          coeff_block[(var1,var2)] = Set{String}()
        end
      end
    end
    size_block = Dict{String, Int64}()
    i_block = 0
    for (block, var_list) in blocks_dict
      i_block += 1
      i_var1 = 1
      size_block[block] = length(var_list)
      jumpX[block] =  @variable(m,[1:size_block[block],1:size_block[block]], base_name= "X_$block", PSD)
      for var1 in var_list
        i_var2 = 1
        for var2 in var_list
          mat_var[(var1,var2)][block] = jumpX[block][i_var1, i_var2]
          coeff_block[(var1,var2)] = union(coeff_block[(var1,var2)],Set{String}([block]))
          i_var2+=1
        end
        i_var1+=1
      end
    end

    xp = Dict{String,JuMP.GenericAffExpr}()
    objctrnames = union(keys(dict_quad_ctr), keys(dict_Bin_ctr), keys(dict_linear_ctr))
    for objctrname in objctrnames
      xp[objctrname] = 0*jumpX["B1"][1,1]
    end

     for (objctrname, dict_ctr) in dict_quad_ctr
       for (vars, value) in dict_ctr
         var1 = vars[1]
         var2 = vars[2]
         nb_block = length(coeff_block[(var1,var2)])
         min_size = Inf
         min_var = first(mat_var[(var1,var2)])[2]
         for (block, var) in mat_var[(var1,var2)]
           if length(blocks_dict[block]) < min_size
             min_size = length(blocks_dict[block])
             min_var = var
           end
         end
         add_to_expression!(xp[objctrname], value*min_var)
       end
     end

     for (objctrname, dict_ctr) in dict_linear_ctr
       for (var, value) in dict_ctr
         add_to_expression!(xp[objctrname], value*λ_and_Sgen_jumpvar[var])
       end
     end

     for (objctrname, value) in dict_constants_ctr
       add_to_expression!(xp[objctrname], value)
     end

    for (objctrname, dict_ctr) in dict_Bin_ctr
      value_binVar =  Dict{String,Float64}()
      Vkk_binVar = Dict{String, GenericAffExpr}()
       for (mono, value) in dict_ctr
         vars = dict_MONO[mono]
         ξ = ""
         Vkk = ""
         name = ""
         abs2_Vkk = 0*jumpX["B1"][1,1]
           if typeof(vars[1]) == Tuple{String,String} || typeof(vars[1]) == Tuple{SubString{String},SubString{String}}
             var1 = vars[1][1]
             var2 = vars[1][2]
             bin_var = vars[2]
           elseif typeof(vars[1]) == String || typeof(vars[1]) == SubString{String}
              bin_var = vars[1]
              var1 = vars[2][1]
              var2 = vars[2][2]
           end
           value_binVar[bin_var] = value
           nb_block = length(coeff_block[(var1,var2)])
           min_size = Inf
           min_var = first(mat_var[(var1,var2)])[2]
           for (block, variable) in mat_var[(var1,var2)]
             if length(blocks_dict[block]) < min_size
               min_size = length(blocks_dict[block])
               min_var = variable
             end
           end
               if !haskey(Vkk_binVar, bin_var)
                 Vkk_binVar[bin_var] = 1* min_var
               else
               add_to_expression!(Vkk_binVar[bin_var], min_var)
              end
        end
        add_to_expression!(xp[objctrname], sum(value*jumpBinvar[bin_var] for (bin_var,value) in value_binVar))
        for (bin_var, abs2_Vkk) in Vkk_binVar
          ξ = jumpBinvar[bin_var]
          name = "ξ_$bin_var"
          constraints_ref[name] = @constraint(m, ξ - abs2_Vkk <= 0)
          uvar[bin_var] = (ξ, abs2_Vkk)
        end
    end

    #add objective
    my_timer = @elapsed @objective(m, Min , xp["OBJ"])

    #add constraints
    for (objctrname, exp) in xp
      if objctrname != "OBJ"
        bounds = dict_bounds_ctr[objctrname]
        if haskey(bounds, "LB")
          lb = bounds["LB"]
        else
          lb = "NONE"
        end
        if haskey(bounds, "UB")
          ub = bounds["UB"]
        else
          ub = "NONE"
        end
        if ub == lb
          constraints_ref[objctrname] = @constraint(m, exp == ub)
        elseif lb != "NONE" && ub != "NONE"
          constraints_ref[objctrname] = @constraint(m, lb <= exp <= ub)
        elseif lb == "NONE" && ub != "NONE"
          constraints_ref[objctrname] = @constraint(m,  exp <= ub)
        elseif lb !="NONE" && ub == "NONE"
          constraints_ref[objctrname] = @constraint(m, lb <= exp)
        end
      end
    end
    #println("initial constraints added to model m")
    #constraints linking common terms in blocks
    nb_coupling_constraints=0
    if NB_BLOCKS > 1
      for i in 1:size(CLIQUE_TREE,1)
        B1 = CLIQUE_TREE[i,1]
        B2 = CLIQUE_TREE[i,2]
        vars_B1 = blocks_dict[B1]
        vars_B2 = blocks_dict[B2]
        common_vars = [ var for var in intersect(vars_B1,vars_B2)]
        for i in 1:length(common_vars)
          var1 = common_vars[i]
          for j in i:length(common_vars)
            var2 = common_vars[j]
            JuMPvar1 = mat_var[(var1,var2)][B1]
            JuMPvar2 = mat_var[(var1,var2)][B2]
            constraints_ref[*("cc$(B1)_$(B2)_$var1","_$var2")]=@constraint(m,JuMPvar1-JuMPvar2==0)
            nb_coupling_constraints+=1
          end
        end
      end
    end
    #print(m)
    optimize!(m)


    println("Objective value : ", JuMP.objective_value(m))
    println("Status : ", JuMP.termination_status(m))

    X_Re, X_Im = construct_approximate_solution(mat_var, blocks_dict, SDP_var_list)

  f = open(joinpath("Mosek_solutions", solution_file), "w")
  write(f,"#Variable    Value \n")
  #
  for (varname, tuple) in uvar
      ξ = tuple[1]
      abs2_Vkk = tuple[2]
      u = JuMP.value(ξ)/JuMP.value(abs2_Vkk)
      #println("$varname : $u")
      write(f, "$varname    $u \n")
  end
  # for (varname, value) in X
  #   write(f, "$varname    $value \n")
  # end
  for (var, val_Re) in X_Re
    value_Re = (JuMP.value(val_Re)) #value_Re = (JuMP.value(X_Re[var]))
    value_Im = (JuMP.value(X_Im[var]))
    write(f, "$(var*"_Re")    $value_Re   \n")
    write(f, "$(var*"_Im")    $value_Im   \n")
  end
  close(f)

  # u_values = [(JuMP.value(tuple[1])/JuMP.value(tuple[2]), varname) for (varname,tuple) in uvar]
  # sort!(u_values, rev=true)
  # f = open(joinpath("fixing4var", solution_file), "w")
  # for i in 1:length(u_values)
  #   tuple = u_values[i]
  #   varname = tuple[2]
  #   if i<=4 && tuple[1] > 0.5
  #     value = 1
  #   else
  #     value = 0
  #   end
  #   write(f, "$varname    $value \n")
  # end
  #
  # close(f)
  return JuMP.objective_value(m), primal_status(m)

end


function solve_SDP(output_instance_path, output_decomposition_path, generation_files_path, FORMULATION, INSTANCE_NAME, generation, flag)
    originalSTDOUT = stdout
    outpath = joinpath("Mosek_runs")
    isdir(outpath) || mkpath(outpath)
    outlog = open(joinpath(outpath,"$(INSTANCE_NAME)_$(generation)_$(flag)_$(FORMULATION).log"), "w")
    redirect_stdout(outlog)
    instance_dat_file_path = output_instance_path
    Pinput_csv_file = joinpath(generation_files_path,"$(INSTANCE_NAME)_$(generation).csv")
    outsolutionpath = joinpath("Mosek_solutions")
    isdir(outsolutionpath) || mkpath(outsolutionpath)
    solution_file = joinpath(outsolutionpath, "$(INSTANCE_NAME)_$(generation)_$(flag).dat")
    λ, Sgen_var_list, SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr,
     dict_MONO, dict_linear_ctr = read_dat_file(instance_dat_file_path)
    cliques_dict, CLIQUE_TREE = read_blocks(output_decomposition_path, FORMULATION, INSTANCE_NAME)
     obj, statut = construct_SDP(cliques_dict, CLIQUE_TREE, Pinput_csv_file, flag, Sgen_var_list, SDP_var_list, Bin_var_list,
          dict_quad_ctr, dict_linear_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, solution_file)
    close(outlog)
    redirect_stdout(originalSTDOUT)
    return obj, statut
end
