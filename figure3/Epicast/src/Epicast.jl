module Epicast

using Statistics, EpicastTables

# ============================================================================ #
abstract type AbstractRunData end
# ============================================================================ #
struct RunData{T} <: AbstractRunData
    demog::FIPSTable{Tract,T,2}
    data::FIPSTable{Tract,T,3}
    fips::Vector{UInt64}
    runno::String
end
rundata(x::RunData, s::AbstractString) = x.data[s]
has_data(x::RunData, s::AbstractString) = EpicastTables.has_var(x.data, s)
demographics(x::RunData, s::AbstractString) = x.demog[s]
n_timepoint(x::RunData) = size(x.data.data, 1)
has_demographic(x::RunData, s::AbstractString) = EpicastTables.has_var(x.demog, s)
function data_groups(x::RunData)
    names = Set{String}()
    for k in keys(x.data.index)
        push!(names, split(k, "_")[1])
    end
    return names
end
column_names(x::RunData) = collect(keys(x.data.var_index))
demographic_names(x::RunData) = collect(keys(x.demog.var_index))
# ============================================================================ #
function run_index(items::AbstractVector{T}) where T
    return Dict{T,Int}(y => x for (x,y) in enumerate(items))
end
# ============================================================================ #
function case_count!(out::AbstractVector, x::RunData, var::AbstractString,
    idx::AbstractVector{<:Integer})

    if has_demographic(x, var)
        # number of new cases relative to total size of that demographic
        out .= new_cases(view(rundata(x, var), idx, :))
        out ./= sum(view(demographics(x, var), idx))
        out .*= 1e5
    else
        out .= total_cases(view(rundata(x, var), idx, :))
    end
    return out
end
# ============================================================================ #
function group_by(x::RunData, name::AbstractString, conv::Integer,
    fcases!::Function=case_count!)

    if has_data(x, name)
        cols = [name]
    else
        cols = filter_vars(x -> startswith(x, name), x.data)
    end

    dat, grp = group_by(x, cols, conv, fcases!)

    return dropdmins(dat, dims=3), grp
end
# ============================================================================ #
function group_by(x::RunData, conv::Integer, fcases!::Function=case_count!)

    cols = sort!(collect(keys(x.data.index)))
    return group_by(x, cols, conv, fcases!)
end
# ============================================================================ #
function group_by(x::RunData, cols::Vector{<:AbstractString}, conv::Integer,
    fcases!::Function=case_count!)

    fips_conv = div.(x.fips, conv)
    unique_grps = sort!(unique(fips_conv))

    out = Array{Float64,3}(undef, n_timepoint(x), length(unique_grps), length(cols))    

    for k in eachindex(unique_grps)
        idx = findall(isequal(unique_grps[k]), fips_conv)
        for j in eachindex(cols)
            fcases!(view(out,:,k,j), x, cols[j], idx)
        end
    end

    return out, unique_grps
end
# ============================================================================ #
function read_runfile_header(io::IO, ::Type{T}=UInt32) where T <: Integer
    nrow = read(io, UInt64)
    ncol = read(io, UInt64)
    n_pt = read(io, UInt64)
    ncol_demog = read(io, UInt64)
    demo_len = read(io, UInt64)
    col_len = read(io, UInt64)

    # @show(Int(nrow), Int(ncol), Int(n_pt), Int(ncol_demog), Int(hdr_len))

    demo_buf = Vector{UInt8}(undef, demo_len)
    read!(io, demo_buf)
    demo_names = map(string, split(String(demo_buf), '\0', keepempty=false))

    col_buf = Vector{UInt8}(undef, col_len)
    read!(io, col_buf)
    col_names = split(String(col_buf), '\0', keepempty=false)

    # FIPS code for each tract
    fips = Vector{UInt64}(undef, nrow)
    read!(io, fips)

    # demographics for of each tract
    demo = Matrix{T}(undef, nrow, ncol_demog)
    read!(io, demo)

    return nrow, ncol, n_pt, col_names, fips,
        FIPSTable{Tract,T,2}(run_index(fips), run_index(demo_names), demo)
end
# ============================================================================ #
function read_runfile(ifile::AbstractString, ::Type{T}=UInt32) where T <: Integer   
    m = match(r".*run_(\d+)\.bin$", ifile)
    m == nothing && error("failed to parse run number, is this a transitions file?")
    run = m[1]
    return open(ifile, "r") do io

        nrow, ncol, n_pt, col_names, fips, demo = read_runfile_header(io, T)

        pos = position(io)
        seekend(io)    
        nbytes = position(io) - pos

        @assert((nbytes % (nrow * ncol * sizeof(UInt32))) == 0,
            "invalid data block size!")
        
        seek(io, pos)

        data = Array{T,3}(undef, nrow, ncol, n_pt)
        read!(io, data)

        return RunData{T}(
            demo,
            FIPSTable{Tract,T,3}(
                run_index(fips),
                run_index(map(string, col_names)),
                permutedims(data, (3,1,2))
            ),
            fips,
            m[1]
        )
    end
end
# ============================================================================ #
struct AgentTransition
    agent_id::UInt64
    location_id::UInt64
    timestep::UInt16
    context::UInt8
    state::UInt8
    variant::UInt8
end
# ---------------------------------------------------------------------------- #
home_state(a::AgentTransition) = a.agent_id >> 58
agent_id(a::AgentTransition) = a.agent_id & ~(UInt64(0b0111111) << 58)
tract_fips(a::AgentTransition) = a.location_id >> 8
tract_community(a::AgentTransition) = a.location_id & 0xff

# (tract, community) in which transition to infected occured
infection_location(a::AgentTransition) = tract_fips(a), tract_community(a)
# ============================================================================ #
function read_agent_transitions(io::IO)
    pos = position(io)
    seekend(io)
    nbyte = position(io) - pos
    seek(io, pos)

    if (nbyte % sizeof(AgentTransition)) != 0
        @warn("File contains inomplete transition packets")
    end

    n_packet = div(nbyte, sizeof(AgentTransition))
    
    data = Vector{AgentTransition}(undef, n_packet)
    read!(io, data)
    
    return data
end
# ============================================================================ #
function RunData(demog::FIPSTable{G,T,2}, data::Vector{AgentTransition},
    fips::Vector{UInt64}, n_pt::Integer, run::AbstractString) where {G,T}

    mp = Dict{UInt64,Int}(id => k for (k,id) in enumerate(fips))

    tmp = zeros(T, n_pt, length(fips), 1)

    # count files do not include index cases
    for d in filter(x -> x.state == 0x01 && x.context != 0xff, data)
        r = mp[tract_fips(d)]
        s = div(d.timestep, 2) + 1

        # count files store data as cumulative sum, so compute that as we go
        view(tmp, s:n_pt, r, 1) .+= 1
    end

    return RunData{T}(
        demog,
        FIPSTable{G,T,3}(run_index(fips), run_index(["total"]), tmp),
        fips,
        run
    )
end
# ============================================================================ #
struct EventData <: AbstractRunData
    demog::FIPSTable{Tract,UInt32,2}
    events::Vector{AgentTransition}
    fips::Vector{UInt64}
    n_pt::UInt64
    run::String
end
# ============================================================================ #
function read_eventfile(::Type{T}, ifile::AbstractString) where T<:AbstractRunData   
    m = match(r".*run_(\d+)\.events\.bin$", ifile)
    m == nothing && error("failed to parse run number, is this a counts file?")
    run = m[1]

    return open(ifile, "r") do io

        nrow, ncol, n_pt, col_names, fips, demog = read_runfile_header(io)

        events = read_agent_transitions(io)

        return T(demog, events, fips, n_pt, m[1])
    end
end

read_eventfile(ifile::AbstractString) = read_eventfile(RunData, ifile)
# ============================================================================ #
total_cases(x::AbstractVector{<:Real}) = x
total_cases(x::AbstractMatrix{<:Real}) = dropdims(sum(x, dims=2),dims=2)
# ============================================================================ #
function new_cases(x::AbstractMatrix{<:Real}, f::Function=sum)
    out = dropdims(f(x, dims=2),dims=2)
    out[2:end] .= diff(out)
    return out
end
function new_cases(x::AbstractVector{<:Real}, f::Function=sum)
    out = copy(x)
    out[2:end] .= diff(out)
    return out
end
mean_new_cases(x) = new_cases(x, mean)
# ============================================================================ #
find_files(dir::AbstractString, re=r".*") = return do_match(dir, re, isfile)
# ---------------------------------------------------------------------------- #
find_directories(dir::AbstractString, re=r".*") = return do_match(dir, re, isdir)
# ============================================================================ #
function do_match(dir::AbstractString, re::Regex, f::Function)
    if !isdir(dir)
        error("Input is not a vaild directory path")
    end
    files = [joinpath(dir, x) for x in readdir(dir)]
    return filter(x->occursin(re, x) && f(x), files)
end
# ============================================================================ #
const STATE_FIPS = Dict(
    1 => "AL", 2 => "AK", 4 => "AZ", 5 => "AR", 6 => "CA", 8 => "CO", 9 => "CT",
    10 => "DE", 11 => "DC", 12 => "FL", 13 => "GA", 15 => "HI", 16 => "ID",
    17 => "IL", 18 => "IN", 19 => "IA", 20 => "KS", 21 => "KY", 22 => "LA",
    23 => "ME", 24 => "MD", 25 => "MA", 26 => "MI", 27 => "MN", 28 => "MS",
    29 => "MO", 30 => "MT", 31 => "NE", 32 => "NV", 33 => "NH", 34 => "NJ",
    35 => "NM", 36 => "NY", 37 => "NC", 38 => "ND", 39 => "OH", 40 => "OK",
    41 => "OR", 42 => "PA", 44 => "RI", 45 => "SC", 46 => "SD", 47 => "TN",
    48 => "TX", 49 => "UT", 50 => "VT", 51 => "VA", 53 => "WA", 54 => "WV",
    55 => "WI", 56 => "WY"
)
# ============================================================================ #
end
