module PltFs
using DrWatson
@quickactivate "OptSitePos"

using Plots
using Dates

export projection_plotter, tunned_plotter, logmessage

DPI=300

function projection_plotter(M, config, weigths, params, flname;show=false,save=true)

  if show == false && save == false
  println(">>>>>>>>>> Skip projection_plotter.")
   return 0
  end

  p = plot(size=(1600, 600), dpi=DPI, leftmargin=10Plots.mm, bottommargin=15Plots.mm,
  tickfontsize=14,xguidefontsize=16,yguidefontsize=16,grid=false)
  for i ∈ eachindex(config)
    if config[i] > 0
      plot!(M[:, 1], ifelse.(M[:, i+1] .≥ params[:clh], i, NaN), color=i,
        linewidth=2.0 * weigths[i, 2], label=false)
      plot!(M[:, 1], ifelse.(params[:cll] .≤  M[:, i+1] .< params[:clh], i, NaN), color=i,
          linewidth= 0.5* weigths[i, 2], label=false)
    end
  end
  ylabel!("Facility Id. number")
  #xticks!([i for i ∈ 0:5:50], [string(i) for i ∈ 0:5:50])
  #yticks!([2*i for i in eachindex(config)],[params[:names][i] for i in eachindex(params[:names])])
  xlabel!("KP [km]")
  if show == true display(p) end
  if save == true savefig(p, plotsdir(savename(flname,params,"pdf"))) end
end



function tunned_plotter(M, config, positions, flname, params::Dict;
  ymin=-105.0, ymax=60.0, hspan=true,show=false,save=true)
  if show == false && save == false
    println(">>>>>>>>>> Skip tunned_plotter.")
    return 0
  end

  p = plot(size=(1600, 600), dpi=DPI, leftmargin=10Plots.mm,
    bottommargin=15Plots.mm, reuse=false,tickfontsize=14,
    xguidefontsize=16,yguidefontsize=16,
    grid=false)
  if hspan
    hspan!([params[:cll], params[:clh]], color=:gray80, label=false)
  end
  for i ∈ eachindex(config)
    if config[i] > 0
      plot!(M[:, 1], M[:, i+1], label=false, color=i,linewidth=0.5)
      if length(positions) == length(config)
        annotate!(positions[i], (40 + ((i % 5) * 5) * (-1)^i), text(params[:names][i], :black, :center, 12))
      end
    end
  end
  #xlims!((-1,27))
  ylims!((ymin, ymax))
  #xticks!([i for i ∈ 0:5:50], [string(i) for i ∈ 0:5:50])
  yticks!([i for i ∈ ymin:15:ymax], [string(i) for i ∈ ymin:15:ymax])
  ylabel!("RX Level [dBm]")
  xlabel!("KP [km]")
  if show == true display(p) end
  if save == true savefig(p, plotsdir(savename(flname,params,"pdf"))) end
end




function logmessage(n)
  # current time
  time = Dates.format(now(UTC), dateformat"yyyy-mm-dd HH:MM:SS")

  # memory the process is using
  memory=round(Sys.maxrss()/1048576, digits=2)
  maxrss = "$(memory) MiB"

  logdata = (;
      n, # iteration n
      maxrss) # lastly the amount of memory being used

  println(">>>>>>>>>> ",savename(time, logdata; connector=" | ", equals=" = ", sort=false, digits=2))
  return memory
end

end
