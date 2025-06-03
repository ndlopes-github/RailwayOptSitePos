#= Copyright (C) 2024
Nuno David Lopes. 
Created:  2024/10/22
Last changed - N. Lopes: 2025/06/03 12:01:24
=#

using DrWatson
@quickactivate "OptSitePos"

using LoggingExtras
using LinearAlgebra
using DataFrames
using CSV
using Gurobi
using JuMP
using SparseArrays

@info "This is the Solvit model for the Antenna Placement Problem"

logger = ConsoleLogger(stderr, Logging.Info; show_limited = false) 
global_logger(logger)
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
	# :LMAXn => 0.01166 * 162.9025,
	:LMAXn => 0.1 * 162.9025,
	# Minimum allowed length for good signal
	:LMINg => 0,
	#:LMINg => 0.85*162.9025,
	#:LMINg => 0.88925 * 162.9025,
	# For restrictions (13)
	# if L=0 do not consider these restrictions
	# in every interval of  length L,
	:L => 0,
	#:L => 5.0,
	# the lengths of the sections without signal do not sum up more than  LMAXnL.
	:LMAXnL => 1.0,


	# Indexes of anchors
	:Ach => [],

	# Data loader/generator and processing options
	# Functions at PreProcessing.jl
	# ["Solvit", LoadJLD2Data"]
	#:preprocess => "Solvit", # inplace preprocessing of data located at data/exp_raw/*.xls 
	:preprocess => "LoadJLD2Data", # Data is preprocessed and should be located at data/exp_pro/solvit_cjMSEnm.jld2
	# Number of antennas
	:nants => 118, #solvit

	#Max number of solutions to search
	:MaxNSol => 20,

	# Images Save and Show
	:ImShow => true,
	:ImSave => false,
)


@info "Parameters $Par"
@info ">>>>>>>>>> Start: Data Processing: "
dataproctime = @elapsed begin
	if Par[:preprocess] == "Solvit"
		# Process specific Solvit data
		# Two files are required for the antennas data
		# One file are required for the priorities and names of antennas
		c, SE, M, nm = pp.Solvit(Par;
			df25 = (datadir("exp_raw", "Douro_coverage_25.xlsx"), "Sheet1"), # Dataframe
			df30 = (datadir("exp_raw", "Douro_coverage_30.xlsx"), "Sheet1"), # Dataframe
			dfw = (datadir("exp_raw", "DouroPriority.xlsx"), "Sheet1"), #Priorities
			# Options for smoothing:
			smoothmethod = ("NoSmoothing",),
			saveprefix = "solvit_")
	elseif Par[:preprocess] == "LoadJLD2Data"
		# Load all data from pre-saved jld2 files.
		# To avoid unnecessary preprocessing of data
		c, SE, M, nm = pp.LoadJLD2Data(Par; filename = "solvit_cjMSEnm.jld2")
	else
		error(">>>>>>>>>> Wrong :preprocess options.")
	end
end
@info ">>>>>>>>>> Data Processing Elapsed time=" dataproctime


pf.projection_plotter(M, ones(Par[:nants]), c, Par, "projection"; show = Par[:ImShow], save = Par[:ImSave])
pf.tunned_plotter(M, ones(Par[:nants]), [], "graphs", Par; ymin = -120, ymax = 45, hspan = true, show = Par[:ImShow], save = Par[:ImSave])


###################################################################################
@info ">>>>>>>>>> Start Optimization Problem construction "
optconsttime = @elapsed begin

	#=
	Partition
	Array SE has the following structure:
	SE = Array{Any}(undef, 0, 4)
	[km -- antena_number --- good/fair(2/2) --- on/off(1/0)]
	=#

	@info ">>>>>>>>>> Partition:"
	SE = SE[sortperm(SE[:, 1]), :]
	counter = copy(SE[:, 4])
	counter .= ifelse.(counter .> 0, 1, -1)
	lngths = zeros(size(SE, 1))
	lngths[1:end-1] .= SE[2:end, 1] .- SE[1:end-1, 1]
	cumsum!(counter, counter)
	SE = hcat(SE, counter)
	SE = hcat(SE, lngths)
	@info "Sample of SE" SE[1:10] 
	println("")
	nosigns = SE[SE[:, 5].==0, 6]
	@info "nosigns lengths array" nosigns
	@info "Sum of nosign lenths" sum(nosigns)
	@info "Maximum of nosigns length" maximum(nosigns)


	# Intervals
	Ip = unique(SE[:, 1])
	# List of Interval lengths
	L = zeros(length(Ip) - 1)
	L[:] .= Ip[2:end] .- Ip[1:end-1]
	global m = length(L)

	@info ">>>>>>>>>> Number of intervals is " m
	maxL = maximum(L[:])
	minL = minimum(L[:])
	@info ">>>>>>>>>> Maximum interval L[i] is" maxL
	@info ">>>>>>>>>> Minimum interval L[i] is" minL

	model = Model(Gurobi.Optimizer)
	hardlimit = 3600 # Max Time For Solver: 1 hour 
	set_optimizer_attribute(model, "TimeLimit", hardlimit)

	# Configure Logging
	set_optimizer_attribute(model, "LogFile", datadir("sims", "gurobi_tune_log.txt")) # Log file path
	set_optimizer_attribute(model, "LogToConsole", 0)
	set_optimizer_attribute(model, "PreSolve", 1)
	set_optimizer_attribute(model, "Heuristics", 0.01)
	set_string_names_on_creation(model, false)
	set_silent(model)

	unregister(model, :x)
	unregister(model, :yg)
	unregister(model, :yf)
	unregister(model, :yn)


	@info "Integer Variables: constraints (11) and (12)"
	@variable(model, 0 ≤ x[1:Par[:nants]] ≤ 1, Int)
	@variable(model, 0 ≤ yg[1:m] ≤ 1, Int)
	@variable(model, 0 ≤ yf[1:m] ≤ 1, Int)
	@variable(model, 0 ≤ yn[1:m] ≤ 1, Int)


	# Objective function (2.1)
	@objective(model, Min, sum((c[j, 2]) * x[j] for j ∈ eachindex(x)))

	# Define and build matrix required for setting the constraints (A)
	function build_A!(A, SE, Ip, fg) #fg=1 "fair" or fg=2 "good"
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

	@info "Model constraint (2)"
	Afg = spzeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
	build_A!(Afg, SE, Ip, 2)
	for i ∈ 1:m, j ∈ 1:Par[:nants]
		if Afg[i, j] == 1.0
			@constraint(model, yg[i] - x[j] ≥ 0)
		end
	end


	@info "Model constraint (3)"
	@constraint(model, [i = 1:m], yg[i] ≤ sum(x[j] * Afg[i, j] for j ∈ 1:Par[:nants]))


	@info "Model constraint (4)"
	Afg = spzeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
	build_A!(Afg, SE, Ip, 1)

	for i ∈ 1:m, j ∈ 1:Par[:nants]
		if Afg[i, j] == 1.0
			@constraint(model, yf[i] - x[j] ≥ 0)
		end
	end

	@info "Model constraint (5)"
	@constraint(model, [i = 1:m], yf[i] ≤ sum(x[j] * Afg[i, j] for j ∈ 1:Par[:nants]))



	@info "Model constraint (6)"
	@constraint(model, [i = 1:m], yg[i] + yf[i] + yn[i] ≥ 1)



	@info "Model constraint (7)"
	@constraint(model, [i=1:m], yg[i] + yn[i] ≤ 1)
	

	@info "Model constraint (8)"
	@constraint(model, [i = 1:m], yf[i] + yn[i] ≤ 1)



	@info "Model constraint (9)"
	@constraint(model, sum(L[i] * yg[i] for i ∈ 1:m) ≥ Par[:LMINg])

	@info "Model constraint (10)"
	@constraint(model, sum(L[i] * yn[i] for i ∈ 1:m) ≤ Par[:LMAXn])

	if Par[:L] > 0
		@info "Model constraint new (13)"
		push!(L, Par[:L] + 1)
		for i ∈ 1:m
			if L[i] > Par[:LMAXnL]
				JuMP.fix(yn[i], 0; force = true) # instead of @constraint(model, yn[i]==0)
			else
				iL = i
				SLK = L[iL]
				while SLK ≤ Par[:L]
					iL += 1
					SLK += L[iL]
				end
				if iL ≤ m
					@constraint(model, sum(L[k] * yn[k] for k ∈ i:(iL-1)) + (Par[:L] - (SLK - L[iL])) * yn[iL] ≤ Par[:LMAXnL])
				else
					@constraint(model, sum(L[k] * yn[k] for k ∈ i:(iL-1)) ≤ Par[:LMAXnL])
				end
			end
		end
		deleteat!(L, length(L))
	end


	@info "Model constraints Anchors"
	for j ∈ Par[:Ach]
		@constraint(model, x[j] == 1.0)
	end

	@info "Model constraints number of antennas"
	if Par[:b] != 0
		@constraint(model, sum(x[j] for j ∈ eachindex(x)) == Par[:b])
	end

end
@info ">>>>>>>>>> Problem construction Elapsed time=" optconsttime


# Allocating space for solutions (0/1)s
sol_xs = zeros(Int, 1, Par[:nants])
sol_ygs = zeros(Int, 1, m)
sol_yfs = zeros(Int, 1, m)
sol_yns = zeros(Int, 1, m)

# reports auxiliar containers
dpts = []
optconsts = []
nintss = []
nantss = []
nvars = []
nconstr = []
stimes = []
obvalues = []
mems = []
tstatus = []

# Solve the optimization model 1st-iteration
@info ">>>>>>>>>> Start Solver: "
t = @elapsed begin
	optimize!(model)
	mem = pf.logmessage(1)

	if termination_status(model) == OPTIMAL
		@info ">>>>>>>>>> Solution is optimal"
	elseif termination_status(model) == TIME_LIMIT && has_values(model)
		@info ">>>>>>>>>> Solution is suboptimal due to a time limit, but a primal solution is available"
	else
		@info ">>>>>>>>>> No solutions"
		error(">>>>>>>>>> The model was not solved correctly.")
	end
	push!(mems, mem)

	nants = Par[:nants]
	push!(nantss, nants)
	push!(nintss, m)

	dpt = dataproctime
	push!(dpts, dpt)
	optconst = optconsttime
	push!(optconsts, optconst)

	stime = solve_time(model)
	push!(stimes, stime)
	@info ">>>>>>>>>> First Solution: solve_time(model) =" stime
	tstat = termination_status(model)
	push!(tstatus, tstat)
	@info ">>>>>>>>>> Termination status = " tstat
	nv = num_variables(model)
	push!(nvars, nv)
	@info ">>>>>>>>>> Number of variables = " nv
	nc = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))
	push!(nconstr, nc)
	@info ">>>>>>>>>> Number of constraints = " nc
	ob = objective_value(model)
	push!(obvalues, ob)
	@info ">>>>>>>>>> Objective value = " ob


	min_cost = round(ob; digits = 2)
	costs = [round(min_cost; digits = 2)] #Rounded objective values

	sol_x = [round(Int, value(x[i])) for i ∈ 1:Par[:nants]]'
	sol_xs[1, :] = sol_x[:]
	sol_yg = [round(Int, value(yg[i])) for i ∈ 1:m]'
	sol_ygs[1, :] = sol_yg[:]
	sol_yf = [round(Int, value(yf[i])) for i ∈ 1:m]'
	sol_yfs[1, :] = sol_yf[:]
	sol_yn = [round(Int, value(yn[i])) for i ∈ 1:m]'
	sol_yns[1, :] = sol_yn[:]


	# Find all the min cost solutions
	global maxnsol = 2
	while result_count(model) ≥ 0 && maxnsol ≤ Par[:MaxNSol]
		@constraint(model, sum(sol_x[j] * x[j] for j ∈ 1:Par[:nants]) ≤ sum(sol_x[:]) - 1)
		optimize!(model)
		global mem = pf.logmessage(maxnsol)

		if termination_status(model) == OPTIMAL
			@info ">>>>>>>>>> Solution is optimal"
		elseif termination_status(model) == TIME_LIMIT && has_values(model)
			@info ">>>>>>>>>> Solution is suboptimal due to a time limit, but a primal solution is available"
		else
			@info ">>>>>>>>>> No more solutions"
			break
		end

		global ob = objective_value(model)
		new_cost = round(ob; digits = 2)
		if new_cost > min_cost
			@info ">>>>>>>>>> Cost Increased, Min Cost = " min_cost "< New Cost =" new_cost
			break
		end
		push!(mems, mem)


		global nants = Par[:nants]
		push!(nantss, nants)
		push!(nintss, m)
		global dpt = dataproctime
		push!(dpts, dpt)
		global optconst = optconsttime
		push!(optconsts, optconst)
		global stime = solve_time(model)
		push!(stimes, stime)
		@info ">>>>>>>>>> " maxnsol "ith Solution: solve_time(model)=" stime
		global tstat = termination_status(model)
		push!(tstatus, tstat)
		@info ">>>>>>>>>> Termination status = " tstat
		global nv = num_variables(model)
		push!(nvars, nv)
		@info ">>>>>>>>>> Number of variables = " nv
		global nc = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))
		push!(nconstr, nc)
		@info ">>>>>>>>>> Number of constraints = " nc
		push!(obvalues, ob)
		@info ">>>>>>>>>> Objective value = " ob

		global sol_x = [round(Int, value(x[i])) for i ∈ 1:Par[:nants]]'
		global sol_xs = [sol_xs; sol_x]
		global sol_yg = [round(Int, value(yg[i])) for i ∈ 1:m]'
		global sol_ygs = [sol_ygs; sol_yg]
		global sol_yf = [round(Int, value(yf[i])) for i ∈ 1:m]'
		global sol_yfs = [sol_yfs; sol_yf]
		global sol_yn = [round(Int, value(yn[i])) for i ∈ 1:m]'
		global sol_yns = [sol_yns; sol_yn]
		global costs = [costs; round(new_cost; digits = 2)]
		global maxnsol += 1
	end
end
@info ">>>>>>>>>> Total Solver Elapsed time=" t

# Plot/show and save pdf files with solutions
for i ∈ 1:size(sol_xs, 1)
	positions = [] # TODO
	pf.tunned_plotter(M, sol_xs[i, :], positions, "graph_sol$(i)", Par; ymax = 55, show = Par[:ImShow], save = Par[:ImSave])
	pf.projection_plotter(M, sol_xs[i, :], c, Par, "projection_sol$(i)"; show = Par[:ImShow], save = Par[:ImSave])
end

@info ">>>>>>>>>> Solutions combinations [x]"
for i ∈ 1:size(sol_xs, 1)
	println(">>>>>>>>>> Combination $(i)")
	display(reshape(sol_xs[i, :], (Int(Par[:nants] / 2), 2)))
end


@info ">>>>>>>>>> Number and names of selected antennas (up to " Par[:nants] " Antennas)"
nonant = sum(sol_xs, dims = 2)
@info nonant
for i ∈ 1:size(sol_xs, 1)
	idx = findall(x -> x == 1, sol_xs[i, :])
	@info ">>>>>>>>>> Solution " i " Antennas list:"
	@info nm[idx]
	@info ">>>>>>>>>> Solution " i " Antennas costs:"
	AC = c[idx, 2]
	@info AC
	@info ">>>>>>>>>> Solution " i " Facilities:"
	@info [(i, count(==(i), AC)) for i in unique(AC)]
end


@info ">>>>>>>>>> Solutions Total Costs (Objective Value) "
@info costs

faircov = (1.0 .- sol_ygs .- sol_yns) * L
@info ">>>>>>>>>> Fair coverage length = " faircov " [km]"
nocov = sol_yns * L
@info ">>>>>>>>>> No coverage lengths = " nocov " [km] . LMAXn = " Par[:LMAXn] "[km]"

@info ">>>>>>>>>> No Coverage extras"
nocovmaxsecs = []
for i ∈ 1:size(sol_yns, 1)
	idx = findall(x -> x == 1, sol_yns[i, :])
	@info ">>>>>>>>>> Solution" i
	considxs = pp.find_consecutive_integers(idx)
	@info ">>>>>>>>>> Sample of consec. indexes " considxs[1:10]
	LengthSums = [sum(L[idxs]) for idxs ∈ considxs[:]]
	@info ">>>>>>>>>> Length of NoCov Sums" length(LengthSums)
	@info ">>>>>>>>>> Sample of Sums" LengthSums[1:10]
	nocovmaxs = maximum(LengthSums)
	@info ">>>>>>>>>> Max of Consecutives yns=" nocovmaxs
	push!(nocovmaxsecs, nocovmaxs)
end


df_sols = DataFrame(sol_xs, nm, makeunique = true)
df_sols[!, "Nants"] .= nantss
df_sols[!, "Nints"] .= nintss
df_sols[!, "NOnAnt"] .= nonant
df_sols[!, "Nvars"] .= nvars
df_sols[!, "Nconstr"] .= nconstr
df_sols[!, "OBvalues"] .= obvalues
df_sols[!, "GoodCov"] .= sol_ygs * L
df_sols[!, "FairCov"] .= faircov
df_sols[!, "NoCov"] .= nocov
df_sols[!, "NoCovMaxL"] .= nocovmaxsecs
df_sols[!, "DPtimes"] .= dpts
df_sols[!, "OptConstT"] .= optconsts
df_sols[!, "SolverTimes"] .= stimes
df_sols[!, "MaxRSS"] .= mems
df_sols[!, "TStatus"] .= tstatus

@info ">>>>>>>>>> Solutions Summary"
@info df_sols[!, Par[:nants]+1:end]


# Write in CSV for easy processing
CSV.write(datadir("sims/solvit", "raw_table_solvit.csv"), df_sols[!, Par[:nants]+1:end-1], delim = ";", append = true)
open(datadir("sims/solvit", "table_solvit.txt"), "a") do file
	write(file, string(df_sols[!, Par[:nants]+1:end]))
	write(file, "\n\n")
end

# Saves a table with all the obtained solutions for one simulation
@info ">>>>>>>>>> Saving Solutions and Reports to CSV"
wsave(datadir("sims/solvit", savename("sol_matrix", Par, "csv")), df_sols)
@info ">>>>>>>>>> END"
