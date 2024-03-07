#=========================================================
  Infeasibility analysis
=========================================================#

const MNT = MutableNamedTuple{(:a1, :q1, :a2, :q2, :d), Tuple{Base.RefValue{Bool}, Base.RefValue{Int}, Base.RefValue{Bool}, Base.RefValue{Int}, Base.RefValue{Float64}}}

# LP model
mutable struct LP
   # variables
   vars::Vector{Int} # id
   active_vars::Vector{Bool}

   # constraints
   cons::Vector{MNT} # a1q1 + a2q2 >= d, true: 1   
   active_cons::Vector{Bool}

   # construtor
   function LP()
      vars = Vector{Int}() 
      active_vars = Vector{Bool}() 
      cons= Vector{MNT}()
      active_cons = Vector{Bool}()
      lp = new(vars, active_vars, cons,active_cons)
      return lp
   end
end

# to name con
@inline con2string(con::MNT) = string(ifelse(con.a1, "q", "-q") , string(con.q1) , ifelse(con.a2, "+q", "-q") , string(con.q2) , ">=" , string(con.d))

# tuple to named tuple
@inline toMNT(a1_::Bool, q1_::Int, a2_::Bool, q2_::Int, d_::Float64) = MutableNamedTuple(a1 = a1_, q1 = q1_, a2 = a2_, q2 = q2_, d = d_)


function printLP(lp)
   print("\n activate vars: ")
   for (i,var) in enumerate(lp.vars)
      if lp.active_vars[i]
         print(var, " ")
      end
   end
   print("\n activate cons: ")
   for (i,con) in enumerate(lp.cons)
      if lp.active_cons[i]
         print("\n ", con2string(con))
      end
   end
end



function merge(var::Int, lcon::MNT, ucon::MNT)
   #print("\n merge: q", var," ", con2string(lcon), " ", con2string(ucon))
   d = lcon.d + ucon.d
   if lcon.q1 == var 
      q1 = lcon.q2
      a1 = lcon.a2
   elseif lcon.q2 == var 
      q1 = lcon.q1
      a1 = lcon.a1
   else 
      assert(false)
   end
   if ucon.q1 == var 
      q2 = ucon.q2
      a2 = ucon.a2
   elseif ucon.q2 == var
      q2 = ucon.q1
      a2 = ucon.a1
   else 
      assert(false)
   end
   return toMNT(a1, q1, a2, q2, d)
end

function sameMNT(con1::MNT, con2::MNT)
   @assert( !(con1.q1 == 0 && con1.a1 == false))
   @assert( !(con1.q2 == 0 && con1.a2 == false))
   @assert( !(con2.q1 == 0 && con2.a1 == false))
   @assert( !(con2.q2 == 0 && con2.a2 == false))
   if con1.q1 == con2.q1 && con1.q2 == con2.q2 && con1.a1 == con2.a1 && con1.a2 == con2.a2
      return true, max(con1.d, con2.d)
   elseif con1.q1 == con2.q2 && con1.q2 == con2.q1 && con1.a1 == con2.a2 && con1.a2 == con2.a1
      return true, max(con1.d, con2.d)
   end
   return false, 0
end

function consofvar(var::Int, lp::LP)
   uppers = Vector{Int}()
   lowers = Vector{Int}()
   for (i, con) in enumerate(lp.cons)
      if lp.active_cons[i] == false
         continue
      end
      if con.q1 == var 
         if con.a1 == true
            push!(lowers, i)
         else
            push!(uppers, i)
         end
      elseif con.q2 == var
         if con.a2 == true
            push!(lowers, i)
         else
            push!(uppers, i)
         end
      end
   end
   return lowers, uppers
end

function checkFeasibleElimination(pi, edge_set, node_set, problem, graph)
   pimap = Dict()
   ef_set = Set()
   for tple in pi
       pimap[tple[1]] = tple[2]
       push!(ef_set, tple[2][1])
   end

   lp = LP()

   var_id = 1
   var_id_map = Dict()
   # construct variables and their bounds
   for ef_id in ef_set
      push!(lp.vars, var_id)
      push!(lp.active_vars, true)
      var_id_map[ef_id] = var_id
      lbd_con = toMNT(true, 0, true, var_id, 0.)
      ubd_con = toMNT(true, 0, false, var_id, -graph.edges[ef_id].length)
      push!(lp.cons, lbd_con)
      push!(lp.cons, ubd_con)
      push!(lp.active_cons, true)
      push!(lp.active_cons, true)
      var_id += 1
   end
   # construct constraints

   for e_id in edge_set
      elen = graph.edges[e_id].length *  (1 + problem.cr_tol) + problem.c_tol
      con = toMNT(true, 0, true, 0, -1.0)
      for i in [:a, :b]
         v_id = graph.edges[e_id].nodes[i]
         var_id = var_id_map[pimap[v_id][1]]
         elen -= problem.dlt - ( problem.d[lor(v_id, graph.edges[pimap[v_id][1]].nodes[pimap[v_id][2]])] +  ifelse(pimap[v_id][2] == :a , 0, (graph.edges[pimap[v_id][1]].length ) )) # big M on x
         if i == :a
            con.a1 = pimap[v_id][2] != :a
            con.q1 = var_id
         else
            con.a2 = pimap[v_id][2] != :a
            con.q2 = var_id
         end
      end
      con.d = elen
      push!(lp.cons, con)
      push!(lp.active_cons, true)
   end


   nactive_vars = length(lp.vars)
   while true
      #printLP(lp)
      # process a single constraint
      redo = false
      for (i, con) in enumerate(lp.cons)
         if lp.active_cons[i] == false
            continue
         end
         # same variables 
         if con.q1 == con.q2
            if con.q1 == 0  # no variables,
               if con.d > problem.c_tol # inconsistent: no variables, 0 > positive
                  return true, con.d 
               else           # deletable
                  lp.active_cons[i] = false
                  redo = true
               end
            elseif lp.active_vars[con.q1] == false # not active
                  continue
            else
               if con.a1 == con.a2 # same sign 
                  con.a1 = true
                  con.q1 = 0
                  con.d /= 2
                  redo = true
               else
                  con.a1 = true
                  con.q1 = 0
                  con.a2 = true
                  con.q2 = 0
                  redo = true
               end
            end
         end
      end
      if redo
         continue
      elseif  nactive_vars == 0
         break
      end
      # proess multiple constraint: tighten bounds and remove redundant constraints
      for (i, con) in enumerate(lp.cons)
         if lp.active_cons[i] == false
            continue
         end
         lencons = length(lp.cons)
         for j in collect(i+1:lencons)
            if lp.active_cons[j] == false
               continue
            end 
            #print("\n", j, "\n", MNT, "\n", typeof(con), "\n", typeof(lp.cons[j]), "\n")
            issame, d = sameMNT(con, lp.cons[j])
            if issame
               lp.active_cons[j] = false
               con.d = d
               redo = true
            end
         end
      end
      if redo
         continue
      end
      # elimination
      for (i, var) in enumerate(lp.vars)
         if lp.active_vars[i] == false
            continue
         end
         lowers, uppers = consofvar(var, lp)
         @assert( !isempty(lowers) )
         @assert( !isempty(uppers) )
         for lcon_id in lowers
            lp.active_cons[lcon_id] = false
            for ucon_id in uppers
               lp.active_cons[ucon_id] = false
               push!(lp.cons, merge(var, lp.cons[lcon_id], lp.cons[ucon_id]))
               push!(lp.active_cons, true)
            end
         end
         lp.active_vars[i] = false
         nactive_vars -= 1
         break
      end
   end

   return false, 0

end


function checkFeasibleLP(model, prev_vars, prev_cons, pi, edge_set, node_set, problem, graph, solver_name, option, time_limit_sec)
   set_silent(model)
   pimap = Dict()
   ef_set = Set()
   for tple in pi
       pimap[tple[1]] = tple[2]
       push!(ef_set, tple[2][1])
   end
   
   # delete 
   for con in prev_cons
       delete(model, con)
       unregister(model, :con)
   end
   for var in prev_vars
       delete(model, var)
       unregister(model, :var)
   end
   if !isempty(prev_vars)
       unregister(model, :q)
       unregister(model, :rv)
   end

   # variables
   @variables(model, begin 
       0 <= q[ef_id in ef_set] <= graph.edges[ef_id].length # edge coordinate variable
       0 <= rv[v_id in node_set] <= problem.Uv[v_id] # residual cover variable
   end) 

   # basic constraints
   @constraints(model, begin 
       [e_id in edge_set], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
       [v_id in node_set], rv[v_id] <=  problem.dlt - ( problem.d[lor(v_id, graph.edges[pimap[v_id][1]].nodes[pimap[v_id][2]])] +  ifelse(pimap[v_id][2] == :a , q[pimap[v_id][1]], (graph.edges[pimap[v_id][1]].length - q[pimap[v_id][1]]) )) # big M on x
   end) 


   # objective
   @objective(model, Min, 0)
   optimize!(model)

   print(model)
   prev_vars = all_variables(model)
   prev_cons = all_constraints(model; include_variable_in_set_constraints = false)

   print(termination_status(model))

   return termination_status(model) == INFEASIBLE, prev_vars, prev_cons
end

function checkFeasibleK1(pi, edge_set, node_set, problem, graph)
   pimap = Dict()
   ef_set = Set()
   for tple in pi
       pimap[tple[1]] = tple[2]
       push!(ef_set, tple[2][1])
   end
   e_id = edge_set[1]
   elen = graph.edges[e_id].length *  (1 - problem.cr_tol) - problem.c_tol

   res = 0
   for v_id in node_set
       mindist = problem.dlt - ( problem.d[lor(v_id, graph.edges[pimap[v_id][1]].nodes[pimap[v_id][2]])])
       res += mindist
       #print(" ", mindist)
   end

   if length(ef_set) == 1
       nodes = collect(node_set)
       v_id1 = nodes[1]
       v_id2 = nodes[2]
       @assert(pimap[v_id1][1] == pimap[v_id2][1])
       issamei = pimap[v_id1][2] == pimap[v_id2][2] 
       if !issamei
           res -= graph.edges[pimap[v_id1][1]].length
       end
   end

   #print("\n", res, " ", elen, " ", length(ef_set), " ", pi, edge_set, node_set, "\n")
   if res < elen
       return true
   else 
       return false
   end

end

function checkAbnormal(pi, edge_set, node_set, problem, graph)
   pimap = Dict()
   ef_set = Set()
   for tple in pi
       pimap[tple[1]] = tple[2]
       push!(ef_set, tple[2][1])
   end
   for e_id in edge_set
       for e_id_ in ef_set
           if e_id in problem.Ec[e_id_]
               return true
           end
       end
   end
   for v_id in node_set
       mindist = problem.dlt - ( problem.d[lor(v_id, graph.edges[pimap[v_id][1]].nodes[pimap[v_id][2]])])
       if mindist < 0
           return true
       end
   end
   return false
end

function enumCoverPatterns(problem::Problem,  option::Option, solver_name::String, K::Int, time_limit_sec)
   graph = problem.prob_graph
   Pi = Vector{Vector{Tuple{Int, Tuple{Int, Symbol}}}}()
   @assert(K <= 3)
   verify = false
   uselp = false
   if uselp || verify
      model = initModel(solver_name, option, time_limit_sec)
      pre_vars = []
      prev_cons = []
   end

   function getNodes(edge_set, prev_node_set)
      total_node_set = Set([[graph.edges[id].nodes[:a] for id in edge_set]; [graph.edges[id].nodes[:b] for id in edge_set]])
      return total_node_set, setdiff(total_node_set, prev_node_set)
   end

   function getPatterns(node_set, prev_patterns)
       node_set = collect(node_set)
       i = 1
       if prev_patterns === nothing
         patterns = [[(node_set[1], enode)] for enode in problem.EIp[node_set[1]]]
         i = 2
       else
         patterns = deepcopy(prev_patterns)
       end
       for vid in node_set[i:end]
           vpatterns = [(vid, enode) for enode in problem.EIp[vid]]
           patterns = [[pattern; tple] for pattern in patterns for tple in vpatterns]
       end
       return patterns
   end

   function isPatternInfeas(pattern, edge_set, node_set)
      if uselp
            infeasible, pre_vars, prev_cons = checkFeasibleLP(model, pre_vars, prev_cons, pattern, edge_set, node_set, problem, graph, solver_name, option, time_limit_sec)
      elseif length(edge_set) == 1
            infeasible = checkFeasibleK1(pattern, edge_set, node_set, problem, graph)
            #infeasible_, pre_vars, prev_cons = checkFeasibleLP(model, pre_vars, prev_cons, pattern, edge_set, node_set, problem, graph, solver_name, option, time_limit_sec)
            #print(infeasible," ", infeasible_)
            #@assert(infeasible == infeasible_)
      else
         if verify 
            infeasible_, pre_vars, prev_cons = checkFeasibleLP(model, pre_vars, prev_cons, pattern, edge_set, node_set, problem, graph, solver_name, option, time_limit_sec)
         end
         infeasible, d = checkFeasibleElimination(pattern, edge_set, node_set, problem, graph)
         if verify 
            print("\n El:", infeasible, " LP:", infeasible_)
            @assert(infeasible == infeasible_ || fabs(d) < problem.c_tol)
         end
      end
      return infeasible
   end

   function isMinPatternInfeas(pattern, edge_set, node_set)
      eset  = (edge_set[2])
      nodes, _ = getNodes(eset, Set())     
      pattern_ = []
      for tple in pattern
         if tple[1] in nodes
            push!(pattern_, tple)
         end
      end
      #print(" ", pattern_, " ", eset, " ", nodes)
      infeas = isPatternInfeas(pattern_, eset, nodes)
      if !infeas && isPatternInfeas(pattern, edge_set, node_set)
         return true
      else
         return false
      end

   end

   
   len_edges = length(graph.edge_ids)
   for i in collect(1:len_edges)
      ei = graph.edge_ids[i]
      edge_set = (ei)
      total_node_seti, _ = getNodes(edge_set, Set())
      patterns = getPatterns(total_node_seti, nothing)
      patterns_i = Vector{Vector{Tuple{Int, Tuple{Int, Symbol}}}}()
      for pattern in patterns
         if checkAbnormal(pattern, edge_set, total_node_seti, problem, graph)
            continue
         end
         infeasible = isPatternInfeas(pattern, edge_set, total_node_seti)
         if infeasible
            push!(Pi, pattern)
         else
            push!(patterns_i, pattern)
         end
      end
      if K == 1
         continue
      end
      for j in collect(i+1:len_edges)
         ej = graph.edge_ids[j]
         if !isAdjacent(graph, ei, ej)
            continue
         end
         edge_set =(ei, ej)
         total_node_setj, diff_node_setj = getNodes(edge_set, total_node_seti)
         patterns = getPatterns(diff_node_setj, patterns_i)
         for pattern in patterns
            if checkAbnormal(pattern, edge_set, total_node_setj, problem, graph)
               continue
            end
            infeasible = isMinPatternInfeas(pattern, edge_set, total_node_setj)
            if infeasible
               push!(Pi, pattern)
            end
         end
      end
   end

   return Pi
end