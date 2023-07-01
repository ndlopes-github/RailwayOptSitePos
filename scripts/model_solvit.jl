using DrWatson
@quickactivate "OptSitePos"

using LinearAlgebra
using DataFrames
using CSV
using HiGHS
using JuMP

include(srcdir("DataPp.jl"))
include(srcdir("Plts.jl"))

import .Pp as pp
import .PltFs as pf


# Parameters
Par = Dict(
  # Model parameters
  # Higher Signal level cut line
  :clh => -80.0,
   # Lower Signal level cut line
  :cll => -95.0,
  # Number of antennas imposed (if 0 constraint is not used);
  :b => 0,
  # Maximum allowed length for no signal
  # 0.01166*162.9025 approx 1.90 KM
  #:LMAXn => 0.01166*162.9025,
  :LMAXn => 0.1*162.9025,
  # Minimum allowed length for good signal
  #:LMINg => 0,
  :LMINg => 0.85*162.9025,
  #:LMINg => 0.88925*162.9025,  
  # For special restrictions (13)
  #if L=0 do not consider these restrictions
  # in every interval of  length L,
  :L => 0.0,
  # the lengths of the sections without signal do not sum up more than  LMAXnL.
  :LMAXnL =>1.0,


  # Indexes of anchors
  :Ach => [ ],

  # Data loader/generator and processing options
  # Functions at PreProcessing.jl
  # ["Solvit", LoadJLD2Data"]
  #:preprocess => "Solvit",
  :preprocess => "LoadJLD2Data",
  # Number of antennas
  :nants=>118, #solvit
  #Max number of solutions to search
  :MaxNSol => 20,

  # Images Save and Show
  :ImShow=>true,
  :ImSave=>true
  )
#@tag!(Par) # DrWatson tag.

println(">>>>>>>>>> Start: Data Processing: ")
dataproctime = @elapsed begin
  # Generic random data is generated
  # Sites with several types of included antennas
  # Work in progress for this one
  # Load specific Solvit data
  # Two files are required for the antennas data
  # One file are required for the priorities and names of antennas
  if Par[:preprocess] == "Solvit"
    c, SE, M, nm = pp.Solvit(Par;
      df25 = (datadir("exp_raw","Douro_coverage_25.xlsx"), "Sheet1"), # Dataframe
      df30 = (datadir("exp_raw","Douro_coverage_30.xlsx"), "Sheet1"), # Dataframe
      dfw = (datadir("exp_raw","DouroPriority.xlsx"), "Sheet1"), #Priorities
      # Options for smoothing:
      smoothmethod = ("NoSmoothing",),
      # smoothmethod = ("ExponentialWeight", 0.1), #Not OK
      # smoothmethod = ("MovingAverage",21) # The same as SavitzkyGolay with degree 1
      # smoothmethod = ("SmoothingSplines", 1.0),
      # smoothmethod =("SavitzkyGolay", (21,5)) #(odd window,polynomial degree)
      saveprefix="raw_")
  # Load all data from pre-saved jld2 files.
  # To avoid unnecessary preprocessing of data
  elseif Par[:preprocess] == "LoadJLD2Data"
      c, SE, M, nm = pp.LoadJLD2Data(Par;loadprefix="raw_")
  else
    error(">>>>>>>>>> Wrong :preprocess options.")
  end
end
println(">>>>>>>>>> Data Processing Elapsed time=", dataproctime)

pf.projection_plotter(M, ones(Par[:nants]), c, Par, "projection";show=Par[:ImShow],save=Par[:ImSave])
pf.tunned_plotter(M,ones(Par[:nants]),[],"graphs",Par; ymin=-120,ymax=45,hspan=true,show=Par[:ImShow],save=Par[:ImSave])


###################################################################################
# Optimization process using JuMP and HiGHS

println(">>>>>>>>>> Start Optimization Problem construction ")
optconsttime = @elapsed begin

  #=
  Partition
  Array SE has the following structure:
  SE = Array{Any}(undef, 0, 4)
  [km -- antena_number --- good/fair(2/2) --- on/off(1/0)]
  =#

  SE = SE[sortperm(SE[:, 1]), :]
  # Intervals
  Ip = unique(SE[:, 1])
  # List of Interval lengths
  L = zeros(length(Ip) - 1)
  L[:] .= Ip[2:end] .- Ip[1:end-1]
  global m=length(L)

  println(">>>>>>>>>> Partition info:")
  println(">>>>>>>>>> Number of intervals is $(m)")
  maxL=maximum(L[:])
  minL=minimum(L[:])
  println(">>>>>>>>>> Maximum interval L[i] is $(maxL)")
  println(">>>>>>>>>> Minimum interval L[i] is $(minL)")


  model = Model(HiGHS.Optimizer)
  set_string_names_on_creation(model, false)
  set_silent(model)
  unregister(model, :x)
  unregister(model, :yg)
  unregister(model, :yf)
  unregister(model, :yn)


  # Integer Variables
  @variable(model, 0 ≤ x[1:Par[:nants]] ≤ 1, Int)
  @variable(model, 0 ≤ yg[1:m] ≤ 1, Int)
  @variable(model, 0 ≤ yf[1:m] ≤ 1, Int)
  @variable(model, 0 ≤ yn[1:m] ≤ 1, Int)


  # Objective function (2.1)
  @objective(model, Min, sum((c[j, 2]) * x[j] for j ∈ eachindex(x)))

  # Define and build matrix required for setting the constraints (A)
  function build_A!(A, SE, I, fg) #fg=1 "fair" or fg=2 "good"
    for (line, i) ∈ enumerate(Ip[1:end-1])
      for (idx, se) ∈ enumerate(SE[:, 1])
        if se == i && Int(SE[idx, 3]) == fg
          A[line, Int(SE[idx, 2])] = SE[idx, 4]
        end
      end
      if line + 1 ≤ m
        A[line+1, :] .= A[line, :]
      end
    end
  end

  # Model constraint (2)
  Ag = zeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
  build_A!(Ag, SE, I, 2)

  for i ∈ 1:m
    for j ∈ 1:Par[:nants]
      if Ag[i, j] == 1.0
        @constraint(model, yg[i] - x[j] ≥ 0)
      end
    end
  end

  # Model constraint (3)
  for i ∈ 1:m
    @constraint(model, yg[i] ≤ sum(x[j] * Ag[i, j] for j ∈ 1:Par[:nants]))
  end

  # Model constraint (4)
  Af = zeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
  build_A!(Af, SE, I, 1)


  for i ∈ 1:m
    for j ∈ 1:Par[:nants]
      if Af[i, j] == 1.0
        @constraint(model, yf[i] - x[j] ≥ 0)
      end
    end
  end

  # Model constraint (5)
  for i ∈ 1:m
    @constraint(model, yf[i] ≤ sum(x[j] * Af[i, j] for j ∈ 1:Par[:nants]))
  end


  # Model constraint (6)
  for i ∈ 1:m
    @constraint(model, yg[i] + yf[i] + yn[i] ≥ 1)
  end


  # Model constraint (7)
  for i ∈ 1:m
    @constraint(model, yg[i] + yn[i] ≤ 1)
  end

  # Model constraint (8)
  for i ∈ 1:m
    @constraint(model, yf[i] + yn[i] ≤ 1)
  end


  # Model constraint (9)
  @constraint(model, sum(L[i] * yg[i] for i ∈ 1:m) ≥ Par[:LMINg])

  # Model constraint (10)
  @constraint(model, sum(L[i] *yn[i] for i ∈ 1:m) ≤ Par[:LMAXn])

  # Model constraint (13)
  if Par[:L]>0
    istar=m
    SL=0
    while SL+L[istar]<Par[:L]
      global istar -= 1
      global SL += L[istar]
    end
    println(">>>>>>>>>> Using restriction (13):  i*=$(istar)")

    i=1
    while i ≤ istar
      local iL=i
      local SL=0
      while SL+L[iL] < Par[:L]
        iL += 1
        SL += L[iL]
      end
      @constraint(model, sum(L[k] *yn[k] for k ∈ i:iL) ≤ Par[:LMAXnL])
      global i += 1
    end
  end


  # Model constraints Anchors
  for j ∈ Par[:Ach]
    @constraint(model, x[j] == 1.0)
  end

  # Model constraints number of antennas
  if Par[:b] != 0
    @constraint(model, sum(x[j] for j ∈ eachindex(x)) == Par[:b])
  end

end
println(">>>>>>>>>> Problem construction Elapsed time=", optconsttime)

# Allocating space for solutions (0/1)s
sol_xs=zeros(Int,1,Par[:nants])
sol_ygs=zeros(Int,1,m)
sol_yfs=zeros(Int,1,m)
sol_yns=zeros(Int,1,m)

# reports auxiliar containers
dpts=[]
optconsts=[]
nintss=[]
nantss=[]
nvars=[]
nconstr=[]
stimes=[]
obvalues=[]
mems=[]
tstatus=[]

# Solve the optimization model 1st-iteration
println(">>>>>>>>>> Start Solver: ")
t = @elapsed begin
  optimize!(model)
  mem=pf.logmessage(1);

  if termination_status(model) == OPTIMAL
    println(">>>>>>>>>> Solution is optimal")
  elseif termination_status(model) == TIME_LIMIT && has_values(model)
    println(">>>>>>>>>> Solution is suboptimal due to a time limit, but a primal solution is available")
  else
    println(">>>>>>>>>> No solutions")
    error(">>>>>>>>>> The model was not solved correctly.")
  end
  push!(mems,mem)

  nants = Par[:nants]; push!(nantss,nants)
  push!(nintss,m)

  dpt = dataproctime; push!(dpts,dpt)
  optconst = optconsttime; push!(optconsts,optconst)

  stime=solve_time(model);push!(stimes,stime)
  println(">>>>>>>>>> First Solution: solve_time(model) =",stime)
  tstat=termination_status(model); push!(tstatus,tstat)
  println(">>>>>>>>>> Termination status = ",tstat)
  nv=num_variables(model);push!(nvars,nv)
  println(">>>>>>>>>> Number of variables = ",nv)
  nc=sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model));push!(nconstr,nc)
  println(">>>>>>>>>> Number of constraints = ",nc)
  ob=objective_value(model);push!(obvalues,ob)
  println(">>>>>>>>>> Objective value = ", ob)


  min_cost = round(ob; digits=2)
  costs = [round(min_cost; digits=2)] #Rounded objective values

  sol_x = [round(Int, value(x[i])) for i ∈ 1:Par[:nants]]'
  sol_xs[1, :] = sol_x[:]
  sol_yg = [round(Int, value(yg[i])) for i ∈ 1:m]'
  sol_ygs[1, :] = sol_yg[:]
  sol_yf = [round(Int, value(yf[i])) for i ∈ 1:m]'
  sol_yfs[1, :] = sol_yf[:]
  sol_yn = [round(Int, value(yn[i])) for i ∈ 1:m]'
  sol_yns[1, :] = sol_yn[:]


  # Find all the min cost solutions
  global maxnsol=2
  while result_count(model) ≥ 0 && maxnsol ≤ Par[:MaxNSol]
    @constraint(model, sum(sol_x[j] * x[j] for j ∈ 1:Par[:nants]) ≤ sum(sol_x[:]) - 1)
    optimize!(model)
    global mem=pf.logmessage(maxnsol);

    if termination_status(model) == OPTIMAL
      println(">>>>>>>>>> Solution is optimal")
    elseif termination_status(model) == TIME_LIMIT && has_values(model)
      println(">>>>>>>>>> Solution is suboptimal due to a time limit, but a primal solution is available")
    else
      println(">>>>>>>>>> No more solutions")
      break
    end

    global ob=objective_value(model)
    new_cost= round(ob; digits=2)
    if new_cost > min_cost
      println(">>>>>>>>>> Cost Increased, Min Cost = ", min_cost, "< New Cost =", new_cost)
      break
    end
    push!(mems,mem)


    # The three following data are constant and are used only for table output coerence
    global nants = Par[:nants]; push!(nantss,nants)
    push!(nintss,m)

    global dpt = dataproctime; push!(dpts,dpt)
    global optconst = optconsttime; push!(optconsts,optconst)
    #####
    global stime=solve_time(model);push!(stimes,stime)
    println(">>>>>>>>>> ", maxnsol, "ith Solution: solve_time(model)=", stime)
    global tstat=termination_status(model);push!(tstatus,tstat)
    println(">>>>>>>>>> Termination status = ",tstat)
    global nv=num_variables(model);push!(nvars,nv)
    println(">>>>>>>>>> Number of variables = ",nv)
    global nc=sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model));push!(nconstr,nc)
    println(">>>>>>>>>> Number of constraints = ",nc)
    push!(obvalues,ob)
    println(">>>>>>>>>> Objective value = ", ob)

    global sol_x = [round(Int, value(x[i])) for i ∈ 1:Par[:nants]]'
    global sol_xs = [sol_xs; sol_x]
    global sol_yg = [round(Int, value(yg[i])) for i ∈ 1:m]'
    global sol_ygs = [sol_ygs; sol_yg]
    global sol_yf = [round(Int, value(yf[i])) for i ∈ 1:m]'
    global sol_yfs = [sol_yfs; sol_yf]
    global sol_yn = [round(Int, value(yn[i])) for i ∈ 1:m]'
    global sol_yns = [sol_yns; sol_yn]
    global costs = [costs; round(new_cost; digits=2)]
    global maxnsol += 1
  end
end
println(">>>>>>>>>> Total Solver Elapsed time=", t)

# Plot/show and save pdf files with solutions
for i ∈ 1:size(sol_xs,1)
    positions=[] # TODO
    pf.tunned_plotter(M,sol_xs[i,:],positions,"graph_sol$(i)",Par;ymax=55,show=Par[:ImShow],save=Par[:ImSave])
    pf.projection_plotter(M, sol_xs[i,:], c, Par, "projection_sol$(i)";show=Par[:ImShow],save=Par[:ImSave])
end

println(">>>>>>>>>> Solutions combinations [x]")
for i ∈ 1:size(sol_xs,1)
  println(">>>>>>>>>> Combination $(i)")
  display(reshape(sol_xs[i,:],(Int(Par[:nants]/2),2)))
end
println(">>>>>>>>>>")

println(">>>>>>>>>> Number and names of selected antennas (up to ", Par[:nants], " Antennas)")
nonant=sum(sol_xs, dims=2)
println(nonant)
for i ∈ 1:size(sol_xs,1)
  idx = findall(x->x==1,sol_xs[i,:])
  println(">>>>>>>>>> Solution ", i, " Antennas list:")
  println(nm[idx])
  println(">>>>>>>>>> Solution ", i, " Antennas costs:")
  AC=c[idx,2]
  println(AC)
  println(">>>>>>>>>> Solution ", i, " Facilities:")
  println([(i, count(==(i), AC)) for i in unique(AC)])
end
println(">>>>>>>>>>")

println(">>>>>>>>>> Solutions Total Costs (Objective Value) ")
println(costs)
println(">>>>>>>>>>")

faircov = (1.0 .- sol_ygs .-sol_yns)*L
println(">>>>>>>>>> Fair coverage length = ",faircov," [km]")
nocov = sol_yns*L
println(">>>>>>>>>> No coverage lengths = ",nocov," [km] . LMAXn = ", Par[:LMAXn], "[km]")

println(">>>>>>>>>> No Coverage extras")
nocovmaxsecs=[]
for i ∈ 1:size(sol_yns,1)
  idx = findall(x->x==1,sol_yns[i,:])
  println(">>>>>>>>>> Solution $(i)")
  considxs=pp.find_consecutive_integers(idx)
  println(">>>>>>>>>> Sample of consec. indexes ",considxs[1:10])
  LengthSums=[sum(L[idxs]) for idxs ∈ considxs[:]]
  nocovmaxs=maximum(LengthSums)
  println(">>>>>>>>>> Max of Consecutives yns=", nocovmaxs)
  push!(nocovmaxsecs,nocovmaxs)
end


df_sols=DataFrame(sol_xs,nm,makeunique=true)
df_sols[!,"Nants"] .= nantss
df_sols[!,"Nints"] .= nintss
df_sols[!,"NOnAnt"] .= nonant
df_sols[!,"Nvars"] .= nvars
df_sols[!,"Nconstr"] .= nconstr
df_sols[!,"OBvalues"] .= obvalues
df_sols[!,"GoodCov"] .= sol_ygs*L
df_sols[!,"FairCov"] .= faircov
df_sols[!,"NoCov"] .= nocov
df_sols[!,"NoCovMaxL"] .= nocovmaxsecs
df_sols[!,"DPtimes"] .= dpts
df_sols[!,"OptConstT"] .= optconsts
stimes[2:end] .= stimes[2: end] .- stimes[1:end-1]
df_sols[!,"Solvertimes"] .= stimes
df_sols[!,"MaxRSS"] .= mems
df_sols[!,"TStatus"] .= tstatus

println(">>>>>>>>>> Solutions Summary")
println(df_sols[!,Par[:nants]+1:end])


# Write in CSV for easy processing
CSV.write(datadir("sims","raw_table.csv"), df_sols[!,Par[:nants]+1:end-1],delim=";", append=true)
open(datadir("sims","table.txt"), "a") do file
  write(file,string(df_sols[!,Par[:nants]+1:end]))
  write(file,"\n\n")
end

# Saves a table with all the obtained solutions for one simulation
println(">>>>>>>>>> Saving Solutions and Reports to CSV")
wsave(datadir("exp_pro",savename("sol_matrix",Par,"csv")),df_sols)
println(">>>>>>>>>>")

println(">>>>>>>>>> END")
