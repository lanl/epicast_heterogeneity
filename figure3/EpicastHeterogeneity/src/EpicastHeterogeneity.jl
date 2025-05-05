module EpicastHeterogeneity

using EpicastGeoplot, PyPlot, Statistics, CSV, Printf
using EpicastGeoplot.Epicast
using EHHelpers
const EH = EHHelpers
# ============================================================================ #
function default_axes!(ax=PyPlot.axes())
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)

    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)

    return ax
end
# ============================================================================ #
function figure3_row(h, ax, data, frame, grp, col)
    
    mx2, time_idc = EpicastGeoplot.add_state_timeseries!(ax[1], data, "total",
        frame, false, "2020-03-20")
    
    ax[1].lines[1].set_color(col)
    
    norm, cm, hp = EpicastGeoplot.add_map!(ax[2], data, "total", frame)

    cb = h.colorbar(
        PyPlot.matplotlib.cm.ScalarMappable(norm=EpicastGeoplot.mpl_norm(norm), cmap=cm),
        cax = ax[3]
    )

    cb.set_label("New cases day $(frame - 1)", fontsize=12)

    default_axes!(ax[4])

    col_lab = column_labels(grp)

    ax[4].bar(col_lab, nanmean(grp.data, dims=1), width=0.9,
        yerr=nanste(grp.data, dims=1), facecolor=col)

    ax[4].set_xticks(col_lab,
        ["white", "black", "asian", "aian", "pi", "other", "multi"], rotation=45)
    ax[4].set_ylabel("Proportion infected", fontsize=14)
    
    ax[2].set_title("")
    ax[1].legend().remove()
end
# ============================================================================ #
function figure3_groupdata(data_file::AbstractString, demo_file::AbstractString,
    t_thr::Integer=typemax(Int))

    data = Epicast.read_eventfile(Epicast.EventData, data_file)

    # ids of agents that transitioned S -> E
    ids = Int64[Epicast.agent_id(x) for x in 
        filter(x -> x.state == 0x01 && x.context != 0xff && x.timestep < t_thr, data.events)]

    # their demographics
    demog = CSV.File(demo_file)

    total_exposed = count(data.events) do x
        x.state == 0x01 && x.context != 0xff
    end

    # count per tract x demographic factor level
    exposed = EH.aggregate_demographic_by(demog[ids .+ 1], :person_race, 1)

    # demographic counts for all agents in each tract in state
    all_agents = EH.aggregate_demographic_by(demog, :person_race, 1)

    out = DictMatrix(exposed.index1, exposed.index2, zeros(size(exposed)))

    for c in column_labels(exposed)
        for r in row_labels(exposed)
            tmp = exposed[r,c] ./ all_agents[r,c]
            out.data[out.index1[r], out.index2[c]] = tmp #isnan(tmp) ? 0.0 : tmp
        end
    end

    return out
end
# ============================================================================ #
function figure3(data_dir::AbstractString)

    data_file = joinpath(data_dir, "figure3", "nm_run_001.events.bin")
    demo_file = joinpath(data_dir, "aux", "nm.csv")
    frame = 93
    data = EpicastGeoplot.geoplot_data(County, data_file)
    grp_data = figure3_groupdata(data_file, demo_file, frame)
    
    h, ax = subplots(3, 4, width_ratios=[0.7, 0.7, 0.05, 0.5])
    h.set_size_inches((9,9.5))

    figure3_row(h, ax[1,:], data, frame, grp_data, "#F8766D")
    h.text(0.5, 0.99, "New Mexico", fontsize=20, ha="center", va="top")

    data_file = joinpath(data_dir, "figure3", "ny_run_001.events.bin")
    demo_file = joinpath(data_dir, "aux", "ny.csv")
    frame = 63
    data = EpicastGeoplot.geoplot_data(County, data_file)
    grp_data = figure3_groupdata(data_file, demo_file, frame)

    figure3_row(h, ax[2,:], data, frame, grp_data, "#00BFC4")
    h.text(0.5, 0.68, "New York", fontsize=20, ha="center", va="top")

    data_file = joinpath(data_dir, "figure3", "fl_run_001.events.bin")
    demo_file = joinpath(data_dir, "aux", "fl.csv")
    frame = 63
    data = EpicastGeoplot.geoplot_data(County, data_file)
    grp_data = figure3_groupdata(data_file, demo_file, frame)

    figure3_row(h, ax[3,:], data, frame, grp_data, "#C77CFF")
    h.text(0.5, 0.36, "Florida", fontsize=20, ha="center", va="top")

    ax[3,1].set_xlabel("Simulation date (month/year)", fontsize=14)

    ax[1,1].set_ylabel("")
    ax[3,1].set_ylabel("")
    ax[1,4].set_ylabel("")
    ax[3,4].set_ylabel("")

    ax[1,1].set_xlabel("")
    ax[2,1].set_xlabel("")

    ax[3,4].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
    ax[3,4].set_ylim(0, 0.015)

    ax[2,4].set_ylim(0, 0.05)

    h.subplots_adjust(top=0.95, right=0.98, left=0.08, hspace=0.5)

    for (max, cax) in zip(ax[:,2],ax[:,3])
        pos = max.get_position()
        wd = (pos.xmax - pos.xmin) * 0.6
        hd = pos.ymax - pos.ymin
        max.set_position([pos.xmin - 0.02, pos.ymin, wd, hd])

        pos2 = cax.get_position()
        wd2 = pos2.xmax - pos2.xmin
        cax.set_position([pos.xmin + wd, pos2.ymin, wd2, hd])
    end

    ax[2,2].set_ylim(ax[2,2].get_ylim() .* [0.99,1.01])

    for (cax, lab) in zip(ax[:,1], ["A.","B.","C."])
        cax.set_ylim(-50, 850)
        pos = cax.get_position()
        h.text(0.075, pos.y1 + 0.005, lab, ha="right", va="bottom", fontsize=24)
    end

    return h, ax
end
# ============================================================================ #
end # module EpicastHeterogeneity
