#= Copyright (C) 2024
Nuno David Lopes.
Created:  2024/10/22
Last changed - N. Lopes:2024/11/22 16:49:02
=#

using DrWatson
@quickactivate "OptSitePos"

using LinearAlgebra
using DataFrames
using CSV
using Gurobi
using JuMP

include(srcdir("DataPp.jl"))
include(srcdir("Plts.jl"))


import .Pp as pp
import .PltFs as pf

#PARAMETERS FOR TEST WITH SIMULATED DataFrame
###############################################
# Use with simsdatagen.sh #
Scale = parse(Int, ARGS[1])
#Scale= 2^5
#Scale=2^5 # Debug
#Scale=2^10
Nsites = Scale * 4
Length = Nsites * 4
Seed = parse(Int, ARGS[2])
#Seed =1
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
  :L => 0.0,
  # the lengths of the sections without signal do not sum up more than  LMAXnL.
  :LMAXnL => 1.0,

  # Indexes of anchors
  :Ach => [],

  # Data loader/generator and processing options
  # Functions at PreProcessing.jl
  :preprocess => "FakeSites",

  # Number of antennas
  :nants => 2 * Nsites, # Number of antennas
  #Max number of solutions to search
  :MaxNSol => 20,

  # Images Save and Show
  :ImShow => false,
  :ImSave => false,

  # Seed for data generation
  :seed => Seed, # 0 for no Seed

  # Initial Max Continous Signal
  :mcs => 15.0 # Lenght for no restriction 
)

println(">>>>>>>>>> Start: Data Processing: ")
dataproctime = @elapsed begin
  # Generic random data is generated
  # Sites with several types of included antennas
  # Work in progress for this one
  if Par[:preprocess] == "FakeSites"
    c, SE, M, nm = pp.FakeSites(Par;
      StartPoint=0.0, # Starting Point in Km
      EndPoint=Length, # Ending Point in Km
      Npoints=10 * Length,
      PrioritiesList=1:3, # Priorities from 1 to 10
      Atypes=[15.0 10.0 25.0; 18.0 12.0 30.0], # Crescent order for all the coordinates
      NFairIs=10, # Max number of fair intervals
      NGoodIs=5, # Max number of good intervals
      save=true # if true: save jld2 preprocessed data
    )
  else
    error(">>>>>>>>>> Wrong :preprocess options.")
  end
end
println(">>>>>>>>>> Data Processing Elapsed time=", dataproctime)


