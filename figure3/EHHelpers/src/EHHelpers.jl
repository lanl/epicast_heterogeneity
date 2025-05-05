module EHHelpers

using CSV, Statistics

export DictMatrix, row_labels, column_labels, nanmean, nanmad, nanstd, nanste,
    nanmedian
# ============================================================================ #
struct DictMatrix{I1,I2,T}
    index1::Dict{I1,Int}
    index2::Dict{I2,Int}
    data::Matrix{T}
end

data(d::DictMatrix) = d.data
Base.size(d::DictMatrix) = Base.size(d.data)
Base.size(d::DictMatrix, k::Integer) = Base.size(d.data, k)
Base.getindex(d::DictMatrix, r, c) = d.data[d.index1[r], d.index2[c]]
Base.getindex(d::DictMatrix, r, ::Colon) = d.data[d.index1[r], 1:size(d.data,2)]
Base.getindex(d::DictMatrix, ::Colon, c) = d.data[1:size(d.data,1), d.index2[c]]
function dict_labels(idx::Dict{T,Int}) where T
    return collect(keys(idx))[sortperm(collect(values(idx)))]
end
row_labels(d::DictMatrix) = dict_labels(d.index1)
column_labels(d::DictMatrix) = dict_labels(d.index2)
# ============================================================================ #
function aggregate_demographic_by(raw::Union{CSV.File,Vector{CSV.Row}},
    demog::Symbol, conv::Integer=1)

    u_tract = sort!(unique(div.(getproperty.(raw, :fips_code), conv)))
    u_demog = sort!(unique(getproperty.(raw, demog)))
    tract_idx = Dict{UInt64,Int}(v => k for (k,v) in enumerate(u_tract))
    demog_idx = Dict{UInt8,Int}(v => k for (k,v) in enumerate(u_demog))

    data = zeros(Int, length(u_tract), length(u_demog))

    for agent in raw
        tract = div(agent.fips_code, conv)
        d = getproperty(agent, demog)
        data[tract_idx[tract], demog_idx[d]] += 1
    end

    return DictMatrix(tract_idx, demog_idx, data)
end
# ============================================================================ #
mutable struct SkipNaN{T<:AbstractVector{<:AbstractFloat}}
    x::T
    _ptr::Int
    _length::Int
end
function skipnan(x::AbstractVector{<:AbstractFloat})
    return SkipNaN(x, 1, length(x) - count(isnan, x))
end
function Base.iterate(s::SkipNaN, k::Integer)
    s._ptr > length(s.x) && return nothing
    ku = findnext(!isnan, s.x, s._ptr)
    ku == nothing && return nothing
    s._ptr = ku + 1
    return s.x[ku], s._ptr
end
Base.iterate(s::SkipNaN) = Base.iterate(s, 1)
Base.IteratorSize(::SkipNaN) = Base.HasLength()
Base.length(s::SkipNaN) = s._length
Base.IteratorEltype(s::SkipNaN) = eltype(s.x)
Base.isdone(s::SkipNaN) = s._ptr > length(s.x)
function Base.collect(s::SkipNaN)
    out = Vector{eltype(s.x)}(undef, s._length)
    k = 1
    @inbounds for j in eachindex(s.x)
        if !isnan(s.x[j])
            setindex!(out, s.x[j], k)
            k += 1
        end
    end
    return out
end
# ============================================================================ #
function nanfun(f::Function, x::AbstractArray{T,2}; dims::Integer=1) where T<:AbstractFloat
    D = abs(dims-3)
    out = Vector{T}(undef,size(x, D))
    k = 1
    for slice in eachslice(x, dims=D)
        out[k] = f(skipnan(slice))
        k += 1
    end
    return out
end
# ============================================================================ #
nanmean(x::AbstractVector{<:AbstractFloat}) = mean(skipnan(x))
nanmean(x::AbstractMatrix{<:AbstractFloat}; dims=1) = nanfun(mean, x, dims=dims)
nanstd(x::AbstractVector{<:AbstractFloat}) = std(skipnan(x))
nanstd(x::AbstractMatrix{<:AbstractFloat}; dims=1) = nanfun(std, x, dims=dims)
@inline ste(itr) = std(itr) / sqrt(length(itr))
nanste(x::AbstractVector{<:AbstractFloat}) = ste(skipnan(x))
nanste(x::AbstractMatrix{<:AbstractFloat}; dims=1) = nanfun(ste, x, dims=dims)
nanmedian(x::AbstractVector{<:AbstractFloat}) = median(skipnan(x))
nanmedian(x::AbstractMatrix{<:AbstractFloat}; dims=1) = nanfun(median, x, dims=dims)
function mad(itr)
    md = median(itr)
    return median((abs.(x .- md) for x in itr))
end
nanmad(x::AbstractVector{<:AbstractFloat}) = mad(skipnan(x))
nanmad(x::AbstractMatrix{<:AbstractFloat}; dims=1) = nanfun(mad, x, dims=dims)
# ============================================================================ #
end # module EHHelpers
