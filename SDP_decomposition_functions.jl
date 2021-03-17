using LightGraphs, MetaGraphs, LinearAlgebra, SparseArrays, SuiteSparse

function load_matpower(filename)
  instance_name = split(filename, '.')[1]

  touch(instance_name*".temp")
  f = open(filename)
  out = open(instance_name*".temp", "w")

  # removing all ';' at end of lines
  while !eof(f)
    line = readline(f)
    if length(line) > 0 && line[1] != '%' && line[1] != 'f'
      s = split(line, ";")
      println(out, s[1])
    end
  end
  close(f)
  close(out)

  data = DelimitedFiles.readdlm(instance_name*".temp")
  rm(instance_name*".temp")
  return data
end

function find_numarray(i_start, data)
  i_debut = i_start
  while !isa(data[i_debut, 1], Int)
    i_debut+=1
  end
  i_fin=i_debut
  while !isa(data[i_fin,1], SubString)
    i_fin += 1
  end
  i_debut, i_fin-1
end

checkfor(data, line_ind, name) = (data[line_ind, 1] == name) || error("Expected ", name, " at line ", line_ind, ", got ", data[line_ind,1], " instead.")

function read_sparsity_pattern(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)

    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    sp = spzeros(nb_bus,nb_bus)

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          sp[orig,dest] = 1
      end
    end

    sp_sym = sp + sp'

    diag = zeros(nb_bus)

    for i=1:nb_bus
        sum_col = sum(sp_sym[i,:])
        diag[i] = Int(sum_col + 1)
    end
    return sp_sym + Diagonal(diag)
    # return sp+sp'+nb_bus*sparse(I, nb_bus, nb_bus)
end

function chordal_ext_cholesky(sparsity_pattern)
    A = sparsity_pattern
    nb_edges_A = (nnz(A) - size(A,1))/2
    #computing cholesky factorisation of A
    F = cholesky(A) #NOTE: order = F.p
    # println(F.p)
    # computing L + LT
    L = sparse(F.L)
    nb_edges_L = nnz(L) - size(A,1)
    nb_added_edges = nb_edges_L - nb_edges_A
    SP = L + L'
    #inverting permutation to get chordal extension of sparsity_pattern
    H = SP[invperm(F.p), invperm(F.p)]
    return H, F.p, nb_added_edges
end

function construct_graph_from_matrix(L)
    n = size(L,1)
    H = MetaGraph(n)
    for i in 1:n
        set_props!(H, i, Dict(:name => "node$i"))
    end
    for i in 1:n
        for j in 1:i
            if L[i,j] != 0
                MetaGraphs.add_edge!(H,i,j)
            end
        end
    end
    return H
end

##########################################################################################################
function cliques_max(G, order)
    n = nv(G)
    C0 = []
    C = C0
    cliques_dict = Dict{Int64,Array{Any,1}}()
    k = 0
    for i=1:n
        v = order[i]
        Nv = neighbors(G,v)
        # println(v)
        C = [Nv[j] for j=1:length(Nv) if Nv[j] ∉ order[1:i]]
        C = [C ; v]
        # println(C)
        if C ⊈ C0
            k +=1
            cliques_dict[k] = C
            C0 = C
        end
    end
    return cliques_dict
end

##############################################################################################################

function weighted_graph(cliques_dict)
    nb_cl = length(cliques_dict)
    G = MetaGraph(nb_cl)
    for i in 1:nb_cl
        for j in 1:i
            inter = length(intersect(cliques_dict[i], cliques_dict[j]))
            if inter > 0
                LightGraphs.add_edge!(G,i,j)
                set_prop!(G, i, j, :weight, inter)
            end
        end
    end
    return G
end

##########################################################################################################

function Prim_algo(G)
    A = MetaGraph(nv(G))
    for i = 1:nv(G)
        set_prop!(A,i, :name, "C$i")
    end
    V = Set([ v for v in LightGraphs.vertices(G)])
    Vprim = Set([1])
    while Vprim != V
        max_weight = 0
        edge_to_add = LightGraphs.Edge(0,0)
        node_to_add = 0
        for edge in LightGraphs.edges(G)
            u = src(edge)
            v = dst(edge)
            if (u ∈ Vprim && v ∉ Vprim)
                weight_uv = get_prop(G,LightGraphs.Edge(u,v), :weight)
                if weight_uv > max_weight
                    max_weight = weight_uv
                    edge_to_add = LightGraphs.Edge(u,v)
                    node_to_add = v
                end
            elseif (u ∉ Vprim && v ∈ Vprim)
                weight_uv = get_prop(G,LightGraphs.Edge(u,v), :weight)
                if weight_uv > max_weight
                    max_weight = weight_uv
                    edge_to_add = LightGraphs.Edge(u,v)
                    node_to_add = u
                end
            end
        end
        LightGraphs.add_edge!(A, edge_to_add)
        set_prop!(A, edge_to_add, :weight, max_weight)
        push!(Vprim, node_to_add)
    end
    return A
end
