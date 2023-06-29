module Pp
using DrWatson
@quickactivate "RailwayOptSitePos"

using XLSX
using DataFrames
using Random
using Distributions

export LoadJLD2Data, Solvit,FakeIntervals

"""
Solvit(Par; df25, df30, dfw, smoothmethod="NotImplemented", save=true, saveprefix="")

This function processes the real data provided by Solvit.

# Arguments
- `Par`: A dictionary containing parameters for the solving operation.
- `df25`: Data frame containing information for antennas of type-1.
- `df30`: Data frame containing information for antennas of type-2.
- `dfw`: Data frame containing the table of priorities.
- `smoothmethod`: (optional-NOT IMPLEMENTED YET) A string representing the smoothing method to be used. Default is "NotImplemented".
- `save`: (optional) A boolean indicating whether to save the results. Default is `true`.
- `saveprefix`: (optional) A string representing the prefix to be used for the saved file. Default is an empty string.

# Output
- `cj`: Cost data frame containing the costs of each antenna.
- `SE`: Array containing information about the partition of the domain.
- `M`: Matrix containing the transferred data from the input data frames.
- `nm`: Array of station names used for distinguishing between 25 and 30 DataFrames.

"""
function Solvit(Par; df25, df30, dfw,
  smoothmethod="NotImplemented",
  save=true,
  saveprefix="")
  # Number of evaluation points in Dataframes
  Npoints = 36420
  ##

  # Loading Table of Priorities
  cj = DataFrame(XLSX.readtable(dfw...))
  cjdict = Dict(cj[i, 1] => parse(Float64, cj[i, 2]) for i ∈ 1:size(cj, 1))

  cjf = 1.2 # Multiplication Factor for costs of antennas30 (Type-2)

  # Data loading and procesing
  df25 = DataFrame(XLSX.readtable(df25...)) #(Type-2)
  df30 = DataFrame(XLSX.readtable(df30...)) #(Type-2)

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
  #println(">>>>>>>>>>> Smoothing method = ", smoothmethod)
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
