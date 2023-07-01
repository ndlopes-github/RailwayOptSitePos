module Pp
using DrWatson
@quickactivate "OptSitePos"

using XLSX
using DataFrames
using Random
using Distributions
using OnlineStats # For Mean and ExponentialWeight
using SmoothingSplines
using SavitzkyGolay

export LoadJLD2Data, Solvit,FakeIntervals


"""
Compute the moving average of an input array using a sliding window of size `m`.

This Docstring was generated using ChatGPT
MovingAverage Function used in the smoothing of real data
source: https://discourse.julialang.org/t/smoothing-noisy-data-using-moving-mean/65329/6

The moving average is computed by sliding a window of size `m` across the input array `A`,
and computing the average of the values within each window. The resulting array has the
same size as the input array, and contains the moving averages.

# Arguments
- `A::AbstractArray`: the input array to compute the moving average for.
- `m::Int`: the size of the sliding window used to compute the moving average. The window
  is centered around each element of the input array, and has a size of `m` elements.

# Output
- `out::AbstractArray`: an array of the same size as the input array, containing the
  computed moving averages.

# Examples
```julia
julia> A = [1, 2, 3, 4, 5]
julia> moving_average(A, 3)
5-element Vector{Float64}:
 1.5
 2.0
 3.0
 4.0
 4.5
"""
function moving_average(A::AbstractArray, m::Int)
  out = similar(A)
  R = CartesianIndices(A)
  Ifirst, Ilast = first(R), last(R)
  I1 = m÷2 * oneunit(Ifirst)
  for I in R
      n, s = 0, zero(eltype(out))
      for J in max(Ifirst, I-I1):min(Ilast, I+I1)
          s += A[J]
          n += 1
      end
      out[I] = s/n
  end
  return out
end


"""
smooth!(M; smoothmethod=("ExponentialWeight", 0.1))

This function takes a matrix M and performs data cleaning and smoothing on each column of the matrix according to the specified method. The function starts by fixing single holes (missing values) in the data. It then applies the selected smoothing method to each column.

Parameters:

    M: a matrix of data to be smoothed
    smoothmethod: a tuple specifying the smoothing method to be used. The first element of the tuple specifies the smoothing method and must be one of the following:
        "NoSmoothing": no smoothing is applied
        "ExponentialWeight": exponential smoothing is applied, where the second element of the tuple specifies the smoothing parameter α. 0 ≤ α ≤ 1, where α=1 means no smoothing.
        "MovingAverage": moving average smoothing is applied, where the second element of the tuple specifies the window size.
        "SmoothingSplines": smoothing splines are applied, where the second element of the tuple specifies the regularization parameter λ.
        "SavitzkyGolay": Savitzky-Golay smoothing is applied, where the second element of the tuple is a tuple specifying the window size and the polynomial degree.

Returns:

    M: the smoothed matrix

Notes:

    The function assumes that missing values in the data are represented by -1111.1111.
    The SavitzkyGolay method works with a window size of 11 and degree of 1.
    Increasing the degree of the SavitzkyGolay method will introduce artificial spikes.

TODO:

    Update the error message to include new smoothing methods.
    SideNote 1: Forecast.jl does not install
    SideNote 2: Loess does not work with NaN
"""
function smooth!(M; smoothmethod=("ExponentialWeight", 0.1))
  # Start by cleaning single holes in data.
  for column in 2:size(M, 2)
    for row in 2:size(M, 1)-1
      # Here we fix one step holes (missing values) in data
      if (
        (M[row-1, column] != -1111.1111) &&
        (M[row+1, column] != -1111.1111) &&
        (M[row, column] == -1111.1111)
      )
        M[row, column] = (M[row-1, column] + M[row+1, column]) / 2.0
      end
    end
    #First and last rows fix.
    if (M[1, column] == -1111.1111) && (M[2, column] != -1111.1111)
      M[1, column] = M[2, column]
    end
    if (M[end, column] == -1111.1111) && (M[end-1, column] != -1111.1111)
      M[end, column] = M[end-1, column]
    end

    if smoothmethod[1] == "ExponentialWeight"
      α = smoothmethod[2] # 0 ≤ α ≤ 1 (no smoothing with α=1)
      o = Mean(weight=ExponentialWeight(α))
      M[:, column] = [value(OnlineStats.fit!(o, yi)) for yi in M[:, column]]
    elseif smoothmethod[1] == "MovingAverage"
      window = smoothmethod[2]
      M[:, column] = moving_average(M[:, column], window)
    elseif smoothmethod[1] == "SmoothingSplines"
      spl = fit(SmoothingSpline, M[:, 1], M[:, column], smoothmethod[2]) # λ=250.0
      M[:, column] = predict(spl) # fitted vector
    elseif smoothmethod[1] == "SavitzkyGolay"
      y = M[:, column]
      sg = savitzky_golay(y, smoothmethod[2][1],smoothmethod[2][2])
      M[:,column] = sg.y
    elseif smoothmethod[1] == "NoSmoothing"
    else
      error(">>>>>>>>>>> Smoothing Method not defined.\n
      Options are: (''NoSmoothing'',) , (''ExponentialWeight'', α), (''MovingAverage'',window))\n
      e.g. α = 0.1 ∈ [0,1] -> Exponential Weight\n
      window=17 -> Moving Average window size.
      Notes:
      NoSmoothing OK
      SavitzkyGolay works with window=11, degree=1
      Increasing degree will introduces artificial spikes

      TODO: UPDATE THIS MESSAGE!
       ")
    end
  end
  M
end



"""
!chatgpt generated: to correct"
Solves an optimization problem based on signal data from two sets of antennas, 25 and 30.
The function loads dataframes containing the signal data, preprocesses them by replacing undefined values,
and checks if they are compatible.
It then extracts a table of priorities for each antenna from an external spreadsheet and reorders it to be compatible with the signal dataframes.
It further renames the stations and combines the signal data from both sets of antennas into a matrix, which is then smoothed using the specified method.
Finally, the function computes a set of subintervals (start and end times) where the signal is above or below certain thresholds for each antenna,
and saves the resulting data (priorities, subintervals, matrix data, and station names) in a file if `save` is true.

# Arguments
- `Par::Dict`: a Dict containing parameters for the optimization problem.
- `df25::String`: a string representing the path to the file containing the signal data from the antennas 25.
- `df30::String`: a string representing the path to the file containing the signal data from the antennas 30.
- `smoothmethod::String = "ExponentialWeight"`: a string representing the method to use for smoothing the signal data. Default is "ExponentialWeight".
- `save::Bool = true`: a boolean representing whether to save the resulting data in a file. Default is true.
- `saveprefix::String = ""`: a string representing a prefix to use for the name of the file where the data will be saved. Default is an empty string.

# Returns
A tuple containing four elements:
- `cj::Matrix{Float64}`: a matrix of priorities for each antenna.
- `SE::Matrix{Any}`: a matrix of subintervals where the signal is above or below certain thresholds for each antenna.
- `M::Matrix{Float64}`: a matrix of the combined signal data from both sets of antennas, smoothed.
- `nm::Vector{String}`: a vector of station names, renamed to distinguish between the antennas 25 and 30.


Function to load and process Solvit data.
Args: Par => Dictionary with the method parameters.
It returns Matrices M (antennas data), Partition SE intervalas, cj priorities and the names (nm) of sites
Solvit Data may be smoothed (smoothing=true) in order to better SE partitions: Work in progress.
This function is taylored to Solvit's Data
e.g.:
c, SE, M, nm = pp.Solvit(Par;
df25 = ("..../data/exp_raw/Douro_coverage_25.xlsx", "Sheet1"), # Dataframe
df30 = ("..../data/exp_raw/Douro_coverage_30.xlsx", "Sheet1"), # Dataframe
dfw = ("..../data/exp_raw/DouroPriority.xlsx", "Sheet1"), #Priorities
smoothing=true)
"""
function Solvit(Par; df25, df30, dfw,
  smoothmethod="ExponentialWeight",
  save=true,
  saveprefix="")
  # Number of evaluation points in Dataframes
  Npoints = 36420
  ##

  # Loading Table of Priorities
  cj = DataFrame(XLSX.readtable(dfw...))
  cjdict = Dict(cj[i, 1] => parse(Float64, cj[i, 2]) for i ∈ 1:size(cj, 1))

  cjf = 1.2 # Multiplication Factor for costs of antennas30

  # Data loading and procesing
  df25 = DataFrame(XLSX.readtable(df25...))
  df30 = DataFrame(XLSX.readtable(df30...))

  # Undefined signal data (no signal) is replaced by -1111.1111
  # In order to simplify partition construction
  df25 = coalesce.(df25, -1111.1111)
  df30 = coalesce.(df30, -1111.1111)

  ## Check if DataFrames are compatible
  if names(df25) != names(df30)
    error(">>>>>>>>>>>>>>>>>>> Incompatible data frames, names of Antennas 25 != names of Antennas 30  <<<<<<<<<<<<<<<<<<<<<<<<<<<  ")
  end

  ## Reorder cj to be compatible with df25 and df30
  ids25 = names(df25, Not("PK"))

  # Priorities are named as in the DataFrames
  cj25 = reshape([ids25; [cjdict[i] for i ∈ ids25]], Par[:nants] ÷ 2, 2)
  # To be modified or generalized
  # For now just multiply by :cjf to the priorities of antennas 25.
  cj30 = reshape([ids25; [cjdict[i] * cjf for i ∈ ids25]], Par[:nants] ÷ 2, 2)
  cj = vcat(cj25, cj30)

  # Display the Costs of each antenna
  display(cj)


  # Rename with shorter ids the stations in order to distinguish
  # 25 and 30 DataFrames station names

  ids25 = [first(name, 5) * "_25" for name ∈ ids25]
  ids30 = names(df30, Not("PK"))
  ids30 = [first(name, 5) * "_30" for name ∈ ids30]

  nm = vcat(ids25, ids30)


  # Transfer data into matrix format
  M25 = Matrix(df25[1:Npoints, :])
  # Transfer data into matrix format
  M30 = Matrix(df30[1:Npoints, :])

  M = hcat(M25, M30[:, 2:end])
  @assert(size(M, 2) - 1 == Par[:nants])

  smooth!(M; smoothmethod=smoothmethod)

  SE = Array{Any}(undef, 0, 4)

  for j ∈ 1:Par[:nants]
    SEs = Array{Any}(undef, 0, 3)
    t = @elapsed begin
      fg = 0
      s_idx = 1
      e_idx = Npoints
      while s_idx < Npoints
        s_idx = findnext(x -> x ≥ Par[:cll], M[:, j+1], s_idx)
        if s_idx === nothing
          break
        end
        fg = ifelse(M[s_idx, j+1] ≥ Par[:clh], 2, 1)
        if fg == 2
          e_idx = findnext(x -> x < Par[:clh], M[:, j+1], s_idx)
          if e_idx === nothing
            e_idx = Npoints
          end
        else
          e_idx = findnext(x -> ((x < Par[:cll]) || (x ≥ Par[:clh])), M[:, j+1], s_idx)
          if e_idx === nothing
            e_idx = Npoints
          end
        end
        SEs = vcat(SEs, [M[s_idx, 1] M[e_idx, 1] fg])
        s_idx = e_idx
      end

      for k ∈ 1:size(SEs, 1)
        SE = vcat(SE, [SEs[k, 1] j SEs[k, 3] 1])
        SE = vcat(SE, [SEs[k, 2] j SEs[k, 3] 0])
      end
    end
    println(">>>>> End Processing Antenna ", j, " Solver Elapsed time=", t)
  end

  if save == true
  wsave(datadir("exp_pro", saveprefix*"cjMSEnm.jld2"),
    Dict("cj" => cj, "M" => M, "SE" => SE, "nm" => nm))
  end
  println(">>>>>>>>>>> Smoothing method = ", smoothmethod)
  cj, SE, M, nm
end



"""
Mconstructor!(M, SEs, Atypes, types, Par, Npoints, k, j)
!Chatgpt generated to correct!

Update a matrix `M` by constructing new entries in the `j`-th column based on
the values of a vector `SEs`, a matrix `Atypes`, a vector `types`, and a
dictionary `Par`.

Arguments:
- `M`: An `Npoints`-by-`j+1` matrix to be updated in-place.
- `SEs`: A `k`-by-`3` matrix with the minimum and maximum values of a range
  and a label (either `1` or `2`) for each range.
- `Atypes`: A `numtypes`-by-`3` matrix of parameters for different types of
  materials.
- `types`: A vector of length `j` containing integer indices into `Atypes`.
- `Par`: A dictionary of parameters used to compute new values for `M`.
- `Npoints`: An integer number of rows in `M`.
- `k`: An integer index into `SEs`.
- `j`: An integer index into the columns of `M` to be updated.

Returns:
- Nothing; `M` is updated in-place.

Examples:
julia> M = rand(3, 3)
julia> SEs = [0.0 1.0 1; 1.0 2.0 2]
julia> Atypes = [1.0 2.0 3.0; 4.0 5.0 6.0]
julia> types = [1, 2, 1]
julia> Par = Dict("cll" => 0.1, "clh" => 0.2)
julia> Npoints = 3
julia> k = 1
julia> j = 2
julia> Mconstructor!(M, SEs, Atypes, types, Par, Npoints, k, j)
"""
function Mconstructor!(M,SEs,Atypes,types,Par,Npoints,k,j)
  for p ∈ 1:Npoints
    if SEs[k, 1] ≤ M[p, 1] ≤ SEs[k, 2]
      if Int(SEs[k, 3]) == 1
        M[p, j+1] = Par[:cll] + 0.75*(Atypes[types[j], 3]/Atypes[end, 3])*(Par[:clh]-Par[:cll])
      elseif Int(SEs[k, 3]) == 2
        M[p, j+1] = Par[:clh] + Atypes[types[j], 3]
      else
        error(" Problemas")
      end
    end
  end
end


"""
DEPRECATED
"""
function FakeIntervals(Par;
  StartPoint=0.0, # Starting Point in Km
  EndPoint=500.0, # Ending Point in Km
  Npoints=50000,
  PrioritiesList=1:10, #Priorities from 1 to 10
  Atypes=[10.0 5.0 25.0; 20.0 7.5 45.0; 35.0 10.0 75.0; 40.0 15.0 100.0], # Fair range Good range Height
  #Atypes = [1.0 0.5 25.0; 2.0 0.75 45.0; 3.5 1.0 75.0; 4.0 1.5 100.0], #Types of intervals
  NFairIs=6,
  NGoodIs=3,
  seed=1123,
  save= true,
  saveprefix="FInt_"
)

  if seed > 0
    Random.seed!(seed)
  end
  Nant = Par[:nants]
  nm = ["A$(i)" for i ∈ 1:Nant]
  Δ = (EndPoint - StartPoint) / (Npoints - 1)


  types = [rand(1:size(Atypes, 1)) for i ∈ 1:Nant]

  #Factor for cj
  Priorities = rand(PrioritiesList, Nant)
  Factors = ones(Nant)

  for i ∈ 1:Nant
    Factors[i] = Atypes[types[i], 3] / (size(Atypes, 1) * minimum(Atypes[:, 3]))
  end

  # Costs are defined by the priorities and by the types of antennas

  cj = hcat(nm, Factors .* Priorities)


  M = (Par[:cll] - 20.0) * ones((Npoints, Nant + 1))
  M[:, 1] = [(StartPoint + (i - 1) * Δ) for i ∈ 1:Npoints]
  SE = Array{Any}(undef, 0, 4) # [km antena number g(2)/f(1) s(1)/e(0)]

  for j ∈ 1:Nant
    SP = StartPoint
    EP = EndPoint
    Nfis = rand(1:NFairIs)
    Ngis = rand(0:NGoodIs)
    Nis = Nfis + Ngis
    SEs = zeros(Nis, 3)
    fgs = shuffle(hcat(ones(Int, Nfis)', 2 * ones(Int, Ngis)'))

    centers = zeros(Nis)
    centers .= StartPoint .+ rand(Uniform(0, 1), Nis) .* (EndPoint - StartPoint)
    sort!(centers)

    for (ifg, fg) ∈ enumerate(fgs)
      s = centers[ifg] - rand(Uniform(0.5, 1))*Atypes[types[j], fg] / 2.0
      e = centers[ifg] + rand(Uniform(0.5, 1))*Atypes[types[j], fg] / 2.0
      if e ≥ EP
        e = EP
      end
      if s ≤ SP
        s = SP
      end
      if s ≥ e
        e=s # This interval will not be considered in the partition
      end

      SEs[ifg, :] = [s e fg]
      SP = e
    end

    for k ∈ 1:Nis
      # Discard degenerate intervals with s >= e
      if SEs[k, 1] ≥ SEs[k, 2]
        continue
      end

      SE = vcat(SE, [SEs[k, 1] j SEs[k, 3] 1])
      SE = vcat(SE, [SEs[k, 2] j SEs[k, 3] 0])

      for p ∈ 1:Npoints
        if SEs[k, 1] ≤ M[p, 1] ≤ SEs[k, 2]
          if Int(SEs[k, 3]) == 1
            M[p, j+1] = (Par[:cll] + Par[:clh]) / 2.0
          elseif Int(SEs[k, 3]) == 2
            M[p, j+1] = Par[:clh] + Atypes[types[j], 3]
          end
        end
      end
    end
  end


  if save == true
  wsave(datadir("exp_pro", saveprefix*"cjMSEnm.jld2"),
        Dict("cj"=> cj,"M"=>M,"SE"=> SE,"nm"=> nm))
  end

  cj, SE, M, nm

end



function FakeSites(Par;
  StartPoint=0.0, # Starting Point in Km
  EndPoint=500.0, # Ending Point in Km
  Npoints=50000,
  PrioritiesList=1:3, #Priorities from 1 to 10
  Atypes=[10.0 6.0 25.0; 12.0 8.0 30.0], # Fair range Good range Height
  NFairIs=6,
  NGoodIs=3,
  seed=1123,
  save=true,
  saveprefix="FSInt_")

  if seed > 0
    Random.seed!(seed)
  end

  Nant = Par[:nants]
  Ntypes = size(Atypes, 1)

  @assert(Ntypes == 2, "Factors for costs defined as in real data, 2 types of Antennas only")
  @assert(Nant % Ntypes == 0, "Number of antennas != number of sites x number of types")

  Nsites = Nant ÷ Ntypes

  # In each site we may install one of each type of antenna
  # Example: [1 1 1 ...2 2 2 ...]
  types = [i for j in 1:Nsites for i ∈ 1:Ntypes]

  #println(types)

  nm = ["A$(j),$(i)" for j in 1:Nsites for i ∈ 1:Ntypes]
  Δ = (EndPoint - StartPoint) / (Npoints - 1)

  #Factor for cj
  Priorities_aux = rand(PrioritiesList, Nsites)
  Priorities = [p for p in Priorities_aux for i in 1:Ntypes]

  #println(Priorities)
  Factors = ones(Nant)

  # The cost factors may be redefined here
  for i ∈ 1:Nant
    # Define as in the Real Data
    # Type 1 => Factor = 1
    # Type 2 => Factor = 1.2
    Factors[i] = (types[i] == 1 ? 1.0 : 1.2 )
  end
  # Costs are defined by the priorities and by the types of antennas
  cj = hcat(nm, Factors .* Priorities)

  if Par[:ImShow] == false && Par[:ImSave] == false
    M = nothing
    println(">>>>>>>>>> Skip M construction (No Plots options)")
  else
    M = (Par[:cll] - 20.0) * ones((Npoints, Nant + 1))
    M[:, 1] = [(StartPoint + (i - 1) * Δ) for i ∈ 1:Npoints]
  end

  SE = Array{Any}(undef, 0, 4) # [km antena number g(2)/f(1) s(1)/e(0)]

  # Define site location for all antennas
  # each site has Ntypes of antennas
  locations = zeros(Nsites)
  locations .= StartPoint .+ rand(Uniform(0, 1), Nsites) .* (EndPoint - StartPoint)

  for tp ∈ 1:Nsites
    Nfis = rand(1:NFairIs)
    Ngis = rand(1:NGoodIs)
    Nis = Nfis + Ngis
    centers = zeros(Nis)


    LocSP = locations[tp] - Atypes[end, 1] / 2.0 #maximum ranges for all types
    LocSP = (LocSP < StartPoint ? StartPoint : LocSP)

    LocEP = locations[tp] + Atypes[end, 1] / 2.0
    LocEP = (LocEP > EndPoint ? EndPoint : LocEP)

    centers .= LocSP .+ rand(Uniform(0, 1), Nis) .* (LocEP - LocSP)
    sort!(centers)

    fgs = shuffle(hcat(ones(Int, Nfis)', 2 * ones(Int, Ngis)'))

    LRndMargin = rand(Uniform(0.7, 1.2), Nis)
    RRndMargin = rand(Uniform(0.7, 1.2), Nis)

    s_s = zeros(Nis, Ntypes, Nsites)
    e_s = zeros(Nis, Ntypes, Nsites)
    SEs = Array{Any}(undef, 0, 3) #zeros(Nis,3)
    ## Larger Reference intervals
    i = Ntypes
    j = tp + (Ntypes - 1) * Nsites
    #println("Antena=", j)
    for (ifg, fg) ∈ enumerate(fgs)
      Ni = (fg == 1 ? Nfis : Ngis)
      s_s[ifg, i, tp] = centers[ifg] - LRndMargin[ifg] * Atypes[end, fg] / (Ni * 2.0)
      e_s[ifg, i, tp] = centers[ifg] + RRndMargin[ifg] * Atypes[end, fg] / (Ni * 2.0)
    end
    SP = StartPoint
    EP = EndPoint

    for (ifg, fg) ∈ enumerate(fgs)
      if e_s[ifg, i, tp] ≥ EP
        e_s[ifg, i, tp] = EP
      end
      if s_s[ifg, i, tp] ≤ SP
        s_s[ifg, i, tp] = SP
      end
      if s_s[ifg, i, tp] ≥ e_s[ifg, i, tp]
        e_s[ifg, i, tp] = s_s[ifg, i, tp]
      end # This interval will not be considered in the partition
      SEs = vcat(SEs, [s_s[ifg, i, tp] e_s[ifg, i, tp] fg])
      #println(SEs[ifg, :])
      SP = e_s[ifg, i, tp]
    end

    #println("Size SEs=", size(SEs, 1))
    for k ∈ 1:size(SEs, 1)
      if M != nothing
        Mconstructor!(M, SEs, Atypes, types, Par, Npoints, k, j)
      end
      if SEs[k, 1] ≥ SEs[k, 2]
        #println("Discard SEs=", SEs[k, :])
        continue
      end # Discard degenerate intervals with s >= e
      SE = vcat(SE, [SEs[k, 1] j SEs[k, 3] 1])
      SE = vcat(SE, [SEs[k, 2] j SEs[k, 3] 0])
    end

    for i ∈ reverse(1:Ntypes-1)
      #SEs=zeros(Nis,3)
      SEs = Array{Any}(undef, 0, 3)
      j = tp + (i - 1) * Nsites
      #println("Antena=", j)
      for (ifg, fg) ∈ enumerate(fgs)
        ###############################################
        # Just need to fine tune the following two lines
        ###############################################
        intlen = (e_s[ifg, i+1, tp] - s_s[ifg, i+1, tp])
        s_s[ifg, i, tp] = s_s[ifg, i+1, tp] + rand(Uniform(0.3, 0.5)) * (1.0 - Atypes[i, fg] / Atypes[end, fg]) * intlen
        e_s[ifg, i, tp] = s_s[ifg, i+1, tp] + rand(Uniform(0.3, 0.5)) * (1.0 + Atypes[i, fg] / Atypes[end, fg]) * intlen

        if s_s[ifg, i, tp] ≥ e_s[ifg, i, tp] # Useless?
          e_s[ifg, i, tp] = s_s[ifg, i, tp] # This interval will not be considered in the partition
        end

        fgaux = 1
        if fg == 2
          fgaux = rand(1:2)
        end

        SEs = vcat(SEs, [s_s[ifg, i, tp] e_s[ifg, i, tp] fgaux])
        #println(SEs[ifg, :])
      end

      #println("Size SEs=", size(SEs, 1))
      for k ∈ 1:size(SEs, 1)
        if M != nothing
          Mconstructor!(M, SEs, Atypes, types, Par, Npoints, k, j)
        end
        if SEs[k, 1] ≥ SEs[k, 2]
          #println("Discard SEs=", SEs[k, :])
          continue
        end # Discard degenerate intervals with s >= e
        SE = vcat(SE, [SEs[k, 1] j SEs[k, 3] 1])
        SE = vcat(SE, [SEs[k, 2] j SEs[k, 3] 0])
      end
    end
  end


  if save == true
    wsave(datadir("exp_pro", saveprefix * "cjMSEnm.jld2"),
      Dict("cj" => cj, "M" => M, "SE" => SE, "nm" => nm))
  end

  cj, SE, M, nm

end


"""
LoadJLD2Data(Par; loadprefix="")

Load data stored in a JLD2 file,
returning the costs coefficients`cj`,
array of partition intervals and info `SE`,
Matrix `M` with kms and values of the signal of each antenna/(set of intervals) for plots,
the list of names/Id's of the antennas/(set of intervals) `nm`.

Arguments:
- `Par`: A dictionary of parameters used to load the data.
- `loadprefix`: Optional string to prepend to the default filename.

Returns:
- `cj`: array of costs coefficients.
- `SE`: array of starts and ends of partition intervals.
- `M`:  Matrix with kms and values of the signal of each antenna/(set of intervals) for plots.
- `nm`: the list of names/Id's of the antennas/(set of intervals).

Examples:
julia> Par = Dict(:clh => -80, :cll => -95, ........)
julia> cj, SE, M, nm = LoadJLD2Data(Par, loadprefix="prefix_")
"""
function LoadJLD2Data(Par; loadprefix="")
  dDict = wload(datadir("exp_pro", loadprefix*"cjMSEnm.jld2"))
  cj = dDict["cj"]
  M = dDict["M"]
  SE = dDict["SE"]
  nm = dDict["nm"]

  
  # Print some data info
  println(">>>>>>>>>> ")
  println("Railway begin=", M[1,1])
  println("Railway end=", M[end,1])
  println("Railway length=", M[end,1]-M[1,1])
  println("Costs informations=", cj)
  println("Base station location and Ids=", nm)
  println(" LoadJLD2Data end")
  println(">>>>>>>>>> ")
  cj, SE, M, nm
end


"""
find_consecutive_integers(arr::AbstractVector{T}) -> Vector{Vector{T}}

Identifies and returns a vector of arrays, each containing a sequence of consecutive integers
found in the input array `arr`.

# Arguments
- `arr::AbstractVector{T}`: An array of integers.

# Returns
- `Vector{Vector{T}}`: A vector of arrays, each containing a sequence of consecutive integers.

# Examples
```julia
julia> arr = [1, 2, 3, 5, 6, 7, 9]
7-element Vector{Int64}:
 1
 2
 3
 5
 6
 7
 9

julia> find_consecutive_integers(arr)
2-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [5, 6, 7]
 This was generated with ChatGPT
"""
function find_consecutive_integers(arr::AbstractVector{T}) where T<:Integer
  consecutive_integers = Vector{T}[]
  n = length(arr)
  i = 1
  while i <= n
      # Find the starting index of a sequence of consecutive integers
      start_index = i
      while i < n && arr[i+1] == arr[i] + 1
          i += 1
      end
      # If a sequence of consecutive integers is found, add it to the result
      # Single indexes are also considered
      if start_index ≤ i
          push!(consecutive_integers, arr[start_index:i])
      end

      i += 1
  end

  return consecutive_integers
end

end
