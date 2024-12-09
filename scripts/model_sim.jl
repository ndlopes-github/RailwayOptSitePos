#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/10/22
Last changed - N. Lopes:2024/12/04 16:39:43
=#

using DrWatson
@quickactivate "OptSitePos"

using Logging
using DataFrames
using CSV
using Dates
using Gurobi
using JuMP
using SparseArrays

@info "This script should be called with loopsims.sh"
@info "Input data should be generated with simsdatagen.sh "

# Set the Input data directory and file prefix
set_dir="sim_set4_2to14_mcs_15_NFair_6_NGood_3/"
file_prefix="cjMSEnm_sim_NFairIs=6_NGoodIs=3_mcs=15.0"

# Set the output log to file
io = open(datadir("sims/paper_table", "logfile_model.txt"), "a")
logger = ConsoleLogger(io, Logging.Info; show_limited=true)
global_logger(logger)

include(srcdir("DataPp.jl"))
include(srcdir("Plts.jl"))


import .Pp as pp
import .PltFs as pf

# INSTANCES 
Instance_I = parse(Int,ARGS[1])
Instance_J = parse(Int,ARGS[2])
NS = [2^(Instance_I+3)] 
IS = [Instance_J] 

for ns ∈ NS, is ∈ IS

    @info "Instance" now() ns is
    instance_dict = @strdict ns is

    Scale = Int(ns * 2^(-3))
    Filename = set_dir*"/"*file_prefix*"_nants=$(ns)_seed=$(is).jld2"
    Nsites = Scale * 4
    Length = Nsites * 4
    ################################################

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
        :LMAXn => 0.1 * Length,
        # Minimum allowed length for good signal
        :LMINg => 0.6 * Length,

        # For special restrictions (13)
        #if L=0 do not consider these restrictions
        # in every interval of  length L,
        :L => 20.0,
        #:L => 0,
        # the lengths of the sections without signal do not sum up more than  LMAXnL.
        :LMAXnL => 15.0,

        # Indexes of anchors
        :Ach => [],

        # Data loader/generator and processing options
        # Functions at PreProcessing.jl
        :preprocess => "LoadJLD2Data",

        # Number of antennas
        :nants => 2 * Nsites, # Number of antennas

        # Images Save and Show
        :ImShow => false,
        :ImSave => false,
        )

    @info Par

    @info ">>>>>>>>>> Start: Data Processing: "

    dataproctime = @elapsed begin
        # Generic random data is generated via sim_instances_generator.jl
        # Sites with several types of included antennas
        if Par[:preprocess] == "LoadJLD2Data"
            c, SE, M, nm = pp.LoadJLD2Data(Par; filename=Filename)
        else
            error(">>>>>>>>>> Wrong :preprocess options.")
        end
    end
    @info ">>>>>>>>>> Data Processing Elapsed time=" dataproctime


    ###################################################################################
    # Optimization process using JuMP and Gurobi

    @info ">>>>>>>>>> Start Optimization Problem construction "

    optconsttime = @elapsed begin
        #=
        Partition
        Array SE has the following structure:
        SE = Array{Any}(undef, 0, 4)
        [km -- antena_number --- good/fair(2/1) --- on/off(1/0)]
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
        println("")
        nosigns = SE[SE[:, 5].==0, 6]
        @info "Sum of nosign lenths" sum(nosigns)
        @info "Maximum of nosigns length" maximum(nosigns)
        @info "Line coverage starts at" SE[1, 1]
        @info "Line coverage end at" SE[end, 1]

        # Intervals
        Ip = unique(SE[:, 1])
        # List of Interval lengths
        L = zeros(length(Ip) - 1)
        L[:] .= Ip[2:end] .- Ip[1:end-1]
        global m = length(L)


        @info ">>>>>>>>>> Number of intervals is" m
        maxL = maximum(L[:])
        minL = minimum(L[:])
        @info ">>>>>>>>>> Maximum interval L[i] is" maxL
        @info ">>>>>>>>>> Minimum interval L[i] is" minL


        model = Model(Gurobi.Optimizer)

        set_string_names_on_creation(model, false)
        #set_silent(model)

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

        # Model constraint (2)
        Ag = spzeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
        build_A!(Ag, SE, Ip, 2)

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
        Af = spzeros(Int, (m, Par[:nants])) # a_ij if antenna j is on in interval i (central antennas)
        build_A!(Af, SE, Ip, 1)


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
        @constraint(model, sum(L[i] * yn[i] for i ∈ 1:m) ≤ Par[:LMAXn])

        # Model constraint (13)
        if Par[:L] > 0
            local k = 1
            local SL = L[k]
            push!(L, 0.0) 
            while k ≤ m
                while SL < Par[:L] && k ≤ m
                    k = k +  1
                    SL = SL + L[k]
                end
                if k ≤ m
                    iL = k
                    SL = L[k]
                    while SL < Par[:L] && k > 1
                        k = k - 1
                        SL = SL + L[k]
                    end
                    i = k
                    @constraint(model, sum(L[k] * yn[k] for k ∈ i:iL) ≤ Par[:LMAXnL])
                    k = k + 1 
                    SL = SL - L[i] + L[k]
                end
            end
            deleteat!(L, length(L)) 
        end
        #

        # Model constraints Anchors
        for j ∈ Par[:Ach]
            @constraint(model, x[j] == 1.0)
        end

        # Model constraints number of antennas
        if Par[:b] != 0
            @constraint(model, sum(x[j] for j ∈ eachindex(x)) == Par[:b])
        end

    end
    @info ">>>>>>>>>> Problem construction Elapsed time=" optconsttime


    # Allocating space for solutions (0/1)s
    global sol_xs = zeros(Int, 1, Par[:nants])
    global sol_ygs = zeros(Int, 1, m)
    global sol_yfs = zeros(Int, 1, m)
    global sol_yns = zeros(Int, 1, m)

    
    dpts = []
    optconsts = []
    nintss = []
    nantss = []
    insts = []
    nvars = []
    nconstr = []
    stimes = []
    obvalues = []
    mems = []
    tstatus = []

    # Start Gurobi Configuration  
    hardlimit = 7200 # Max Time For Solver: 120 mnts 
    set_optimizer_attribute(model, "TimeLimit", hardlimit)
    # Configure Logging
    gurobi_log_instance = savename("gurobi_log", instance_dict, "txt")
    set_optimizer_attribute(model, "LogFile", datadir("sims/paper_table", gurobi_log_instance)) # Log file path
    set_optimizer_attribute(model, "OutputFlag", 1) # Output info to log
    set_optimizer_attribute(model, "LogToConsole", 1) # 0: Disable console logging, 1: Enable
    set_optimizer_attribute(model, "Heuristics", 0.01) # Ranges from 0 to 1, where 0 disables heuristics and 1 maximizes heuristic effort. The default is 0.05 (5% of effort).
    set_optimizer_attribute(model, "PreSolve", 2) # Controls the presolve level. (range -1 to 3; -1 = auto, 0 = no presolve, 1 = conservative, 2 = aggressive)
    # End Gurobi Configuration


    # Solve the optimization model 1st-iteration
    @info ">>>>>>>>>> Start Solver: "
    t = @elapsed begin

        optimize!(model)
        
        if !is_solved_and_feasible(model; dual=false)
            @warn(
                """
                The model was not solved correctly:
                termination_status : $(termination_status(model))
                primal_status      : $(primal_status(model))
                dual_status        : $(dual_status(model))
                raw_status         : $(raw_status(model))
                """,
            )
            # If infeasible or not solved continue to next instance
            continue 
        end

        mem = pf.logmessage(1)
        push!(mems, mem)

        
        global nants = Par[:nants]
        push!(nantss, nants)
        push!(insts,is)
        push!(nintss, m)
        
        global dpt = dataproctime
        push!(dpts, dpt)
        global optconst = optconsttime
        push!(optconsts, optconst)

        global stime = solve_time(model)
        push!(stimes, stime)
        @info ">>>>>>>>>> First Solution: solve_time(model) $(stime)"
        global tstat = termination_status(model)
        push!(tstatus, tstat)
        @info ">>>>>>>>>> Termination status = " tstat
        global nv = num_variables(model)
        push!(nvars, nv)
        @info ">>>>>>>>>> Number of variables = " nv
        global nc = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))
        push!(nconstr, nc)
        @info ">>>>>>>>>> Number of constraints = " nc
        global ob = objective_value(model)
        push!(obvalues, ob)
        @info ">>>>>>>>>> Objective value = " ob


        min_cost = round(ob; digits=2)
        global costs = [round(min_cost; digits=2)] #Rounded objective values
        global sol_x = [round(Int, value(x[i])) for i ∈ 1:Par[:nants]]'
        sol_xs[1, :] = sol_x[:]
        global sol_yg = [round(Int, value(yg[i])) for i ∈ 1:m]'
        sol_ygs[1, :] = sol_yg[:]     
        global sol_yf = [round(Int, value(yf[i])) for i ∈ 1:m]'
        sol_yfs[1, :] = sol_yf[:]
        global sol_yn = [round(Int, value(yn[i])) for i ∈ 1:m]'
        sol_yns[1, :] = sol_yn[:]

    end
    @info ">>>>>>>>>> Total Solver Elapsed time=" t

    @info "solution summary" solution_summary(model; result=1, verbose=true)

    @info ">>>>>>>>>> Number and names of selected antennas (up to " Par[:nants] " Antennas)"
    nonant = sum(sol_xs, dims=2)
    @info nonant
    for i ∈ 1:size(sol_xs, 1)
        idx = findall(x -> x == 1, sol_xs[i, :])
        AC = c[idx, 2]
        @info [(i, count(==(i), AC)) for i in unique(AC)]
    end

    @info ">>>>>>>>>> Solutions Total Costs (Objective Value) "
    @info costs

    faircov = (1.0 .- sol_ygs .- sol_yns) * L
    @info ">>>>>>>>>> Fair coverage length = " faircov " [km]"
    nocov = sol_yns * L
    @info ">>>>>>>>>> No coverage lengths = $(nocov)  [km] ."
    @info ">>>>>>>>>> No Coverage extras"

    nocovmaxsecs = []
    for i ∈ 1:size(sol_yns, 1)
        idx = findall(x -> x == 1, sol_yns[i, :])
        @info ">>>>>>>>>> Solution " i
        considxs = pp.find_consecutive_integers(idx)
        if length(considxs) > 0
            # @info ">>>>>>>>>> Sample of consec. indexes " considxs[1:end]
            LengthSums = [sum(L[idxs]) for idxs ∈ considxs[:]]
            nocovmaxs = maximum(LengthSums)
            @info ">>>>>>>>>> Max of Consecutives yns =  $(nocovmaxs)"
            push!(nocovmaxsecs, nocovmaxs)
        end
    end


    df_sols = DataFrame(sol_xs, nm, makeunique=true)
    df_sols[!, "Nants"] .= nantss
    df_sols[!, "Inst"] .= insts
    df_sols[!, "Nints"] .= nintss
    df_sols[!, "NOnAnt"] .= nonant
    df_sols[!, "Nvars"] .= nvars
    df_sols[!, "Nconstr"] .= nconstr
    df_sols[!, "OBvalues"] .= obvalues
    df_sols[!, "GoodCov"] .= sol_ygs * L
    df_sols[!, "FairCov"] .= faircov
    df_sols[!, "NoCov"] .= nocov

    if length(nocovmaxsecs) > 0
        df_sols[!, "NoCovMaxL"] .= nocovmaxsecs
    else
        nocovmaxsecs = zeros(length(nocov))
        df_sols[!, "NoCovMaxL"] .= nocovmaxsecs
    end

    df_sols[!, "DPtimes"] .= dpts
    df_sols[!, "OptConstT"] .= optconsts
    df_sols[!, "SolverTimes"] .= stimes
    df_sols[!, "MaxRSS"] .= mems
    df_sols[!, "TStatus"] .= tstatus

    @info ">>>>>>>>>> Solutions Summary"
    @info df_sols[!, Par[:nants]+1:end]


    # Write in CSV for easy processing
    CSV.write(datadir("sims/paper_table", "raw_table.csv"), df_sols[!, Par[:nants]+1:end-1], delim=";", append=true)
    open(datadir("sims/paper_table", "table.txt"), "a") do file
        write(file, string(df_sols[!, Par[:nants]+1:end]))
        write(file, "\n\n")
    end

    @info ">>>>>>>>>> END"
    flush(io)
end
close(io)
