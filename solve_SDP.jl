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

function read_blocks(DATA_PATH, FORMULATION, INSTANCE_NAME)
    if FORMULATION == "NO_BLOCKS"
        BLOCKS=readdlm(joinpath(DATA_PATH,"blocks_cholesky", INSTANCE_NAME*"_sdp_blocks.txt"))
        VAR = BLOCKS[:,2]
        B1 = ["B1" for i=1:size(BLOCKS,1)]
        BLOCKS = [B1 VAR]
        CLIQUE_TREE = []
    else
        BLOCKS=readdlm(joinpath(DATA_PATH, "blocks_$FORMULATION",INSTANCE_NAME*"_sdp_blocks.txt"))
        CLIQUE_TREE = readdlm(joinpath(DATA_PATH, "cliquetree_$FORMULATION",INSTANCE_NAME*"_sdp_cliquetree.txt"))
    end
    cliques_dict = Dict{String,Array{String}}()
    for i in 1:size(BLOCKS,1)
        block = BLOCKS[i,1]
        var = BLOCKS[i,2]
            if !haskey(cliques_dict, block)
                cliques_dict[block] = [var]
            else
                push!(cliques_dict[block], var)
            end
    end
    return cliques_dict, CLIQUE_TREE
end

function construct_approximate_solution(mat_var, blocks_dict, SDP_var_list)
  #complex blocks
  # f = open("sol_C.csv", "w")
  eigenvectors_max = Dict{String,Array{Complex{Float64}}}()
  dict_C_blocks = Dict{String, Array{Complex{Float64}}}()
  blocks_dict_C = Dict{String, Array{String}}()
  for (block, list_var) in blocks_dict
    nb_var_R = length(list_var)
    nb_var_C = Int(nb_var_R/2)
    list_var_C = [var[1:end-3] for var in list_var if var[end-1:end]=="Re"]
    # println(list_var_C)
    W = Array{Complex{Float64}}(undef, nb_var_C, nb_var_C)
    for k in 1:nb_var_C
      var1_C =  list_var_C[k]
      var1_Re = var1_C*"_Re"
      var1_Im = var1_C*"_Im"
      for l in k:nb_var_C
        var2_C =  list_var_C[l]
        var2_Re = var2_C*"_Re"
        var2_Im = var2_C*"_Im"
        if k==l
          W[k,l] = JuMP.value(mat_var[(var1_Re,var2_Re)][block]) + JuMP.value(mat_var[(var1_Im,var2_Im)][block])
        else
          W[k,l] = JuMP.value(mat_var[(var1_Re,var2_Re)][block]) + JuMP.value(mat_var[(var1_Im,var2_Im)][block])+ im*(JuMP.value(mat_var[(var1_Im,var2_Re)][block])-JuMP.value(mat_var[(var1_Re,var2_Im)][block]))
          W[l,k] = JuMP.value(mat_var[(var1_Re,var2_Re)][block]) +  JuMP.value(mat_var[(var1_Im,var2_Im)][block])- im*(JuMP.value(mat_var[(var1_Im,var2_Re)][block])-JuMP.value(mat_var[(var1_Re,var2_Im)][block]))
        end
      end
    end
    dict_C_blocks[block] = W
    # println(ishermitian(W))
    # println(dict_C_blocks[block])
    X_block = dict_C_blocks[block]
    K = size(X_block, 1)
    ef = eigen(X_block)
    ef.values[K:K]       #largest k
    eigenvectors_max[block] = ef.vectors[:, K:K]
    # println(eigenvectors_max)
    # norm_eig = norm(ef.vectors[:, K:K])
    # println(norm_eig)
    for i=1:K
      var = list_var_C[i]
      value = eigenvectors_max[block][i]
      # write(f, "$block ; $var ; $value  \n")
    end
    blocks_dict_C[block] = list_var_C
  end
  list_blocks = [block for block in keys(blocks_dict_C)]
  NB_BLOCKS = length(list_blocks)
  sync_pb = Model(with_optimizer(Xpress.Optimizer))
  @variable(sync_pb, θ[1:NB_BLOCKS], lower_bound=0, upper_bound=2*pi)
  obj = @expression(sync_pb, 0.0*θ[1]^2)
  for i in 1:NB_BLOCKS
    block1 = list_blocks[i]
    list_var_1 = blocks_dict_C[block1]
    p1 = eigenvectors_max[block1]
    for j in (i+1):NB_BLOCKS
      block2 = list_blocks[j]
      list_var_2 = blocks_dict_C[block2]
      p2 = eigenvectors_max[block2]
      for k in 1:length(list_var_1)
        var1 = list_var_1[k]
        for l in 1:length(list_var_2)
          var2 = list_var_2[l]
          if var1 == var2
            add_to_expression!(obj, (angle(p1[k])+θ[i]-angle(p2[l])-θ[j])^2)
          end
        end
      end
    end
  end
  @objective(sync_pb, Min, obj)
  optimize!(sync_pb)
  findX = Model(with_optimizer(Xpress.Optimizer))
  X_Re = Dict{String, VariableRef}()
  X_Im = Dict{String, VariableRef}()
  list_var_C = [var[1:end-3] for var in SDP_var_list if var[end-1:end]=="Re"]
  # println(list_var_C)
  for var in list_var_C
    X_Re[var] = @variable(findX, base_name="X_Re_$var")
    X_Im[var] = @variable(findX, base_name="X_Im_$var")
  end

  obj = @expression(findX, 0.0*X_Re["VOLT_1"]^2)
  # f = open("sol_sync_pb.csv", "w")
  for i in 1:NB_BLOCKS
    block = list_blocks[i]
    θ_block = JuMP.value(θ[i])
    p = exp(im*θ_block).*eigenvectors_max[block]
    i = 1
    for var in blocks_dict_C[block]
      value = p[i]
      # write(f, "$block ; $var ; $(abs(eigenvectors_max[block][i])) ; $(angle(eigenvectors_max[block][i])) ; $(abs(value)) ; $(angle(value)) \n")
      add_to_expression!(obj, (X_Re[var]-real(value))^2+(X_Im[var]-imag(value))^2)

      i +=1
    end
  end
  @objective(findX, Min, obj)
  optimize!(findX)
  for var in list_var_C
    value = (JuMP.value(X_Re[var])) + im*(JuMP.value(X_Im[var]))
  end
  return X_Re, X_Im
end

function construct_SDP(blocks_dict, CLIQUE_TREE, Pinput_csv_file, flag, Sgen_var_list, SDP_var_list, Bin_var_list,
   dict_quad_ctr, dict_linear_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, solution_file)

   lines = readlines(Pinput_csv_file)
   NB_BLOCKS = length(blocks_dict)
    #initialize
    if solution_file == "" #SDP in B&B
        m = Model(with_optimizer(Mosek.Optimizer, LOG=0))
    else
        m = Model(with_optimizer(Mosek.Optimizer))
    end    
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

    value_bins = Dict{String, Float64}()
    for (varname, tuple) in uvar
        ξ = tuple[1]
        abs2_Vkk = tuple[2]
        u = JuMP.value(ξ)/JuMP.value(abs2_Vkk)
        value_bins[varname] = u
    end
    value_SDP_var = Dict{String, Array{Float64,2}}()
    for (block, X) in jumpX
      value_SDP_var[block] = JuMP.value.(X)
    end

    if solution_file != ""
        X_Re, X_Im = construct_approximate_solution(mat_var, blocks_dict, SDP_var_list)
        isdir("Mosek_solutions") || mkpath("Mosek_solutions")
          f = open(joinpath("Mosek_solutions", solution_file), "w")
          write(f,"#Variable    Value \n")
          for (varname, tuple) in uvar
              ξ = tuple[1]
              abs2_Vkk = tuple[2]
              u = JuMP.value(ξ)/JuMP.value(abs2_Vkk)
              write(f, "$varname    $u \n")
          end
          for (var, val_Re) in X_Re
            value_Re = (JuMP.value(val_Re)) #value_Re = (JuMP.value(X_Re[var]))
            value_Im = (JuMP.value(X_Im[var]))
            write(f, "$(var*"_Re")    $value_Re   \n")
            write(f, "$(var*"_Im")    $value_Im   \n")
          end
          close(f)
      end

  return value_bins, value_SDP_var, JuMP.objective_value(m), primal_status(m)

end



function solve_SDP(ROPF, flag)
    output_instance_path = ROPF.output_instance_path
    output_decomposition_path = ROPF.output_decomposition_path
    generation_files_path = ROPF.generation_files_path
    FORMULATION = ROPF.decomposition
    INSTANCE_NAME = ROPF.instance_name
    generation = ROPF.generation
    originalSTDOUT = stdout
    outpath = joinpath("Mosek_runs")
    isdir(outpath) || mkpath(outpath)
    outlog = open(joinpath(outpath,"$(INSTANCE_NAME)_$(generation)_$(flag)_$(FORMULATION).log"), "w")
    redirect_stdout(outlog)
    instance_dat_file_path = joinpath(output_instance_path, "$(INSTANCE_NAME).dat")
    Pinput_csv_file = joinpath(generation_files_path,"$(INSTANCE_NAME)_$(generation).csv")
    solution_file = joinpath("$(INSTANCE_NAME)_$(generation)_$(flag).dat")
    λ, Sgen_var_list, SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr,
     dict_MONO, dict_linear_ctr = read_dat_file(instance_dat_file_path)
    cliques_dict, CLIQUE_TREE = read_blocks(output_decomposition_path, FORMULATION, INSTANCE_NAME)
     value_bins, value_SDP_var, obj, statut = construct_SDP(cliques_dict, CLIQUE_TREE, Pinput_csv_file, flag, Sgen_var_list, SDP_var_list, Bin_var_list,
          dict_quad_ctr, dict_linear_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, solution_file)
    close(outlog)
    redirect_stdout(originalSTDOUT)
    return obj, statut
end
