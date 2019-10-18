struct RiskPriors{M <: ILM}
  sparks::Union{Nothing, Vector{UnivariateDistribution}}
  susceptibility::Union{Nothing, Vector{UnivariateDistribution}}
  infectivity::Union{Nothing, Vector{UnivariateDistribution}}
  transmissibility::Union{Nothing, Vector{UnivariateDistribution}}
  latency::Union{Nothing, Vector{UnivariateDistribution}}
  removal::Union{Nothing, Vector{UnivariateDistribution}}

  function RiskPriors{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F, Ωr::F) where {
                         F <: Vector{UnivariateDistribution}, S <: SEIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, Ωl)
  end

  function RiskPriors{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F) where {
                         F <: Vector{UnivariateDistribution}, S <: SEI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, nothing)
  end

  function RiskPriors{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωr::F) where {
                         F <: Vector{UnivariateDistribution}, S <: SIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, Ωr)
  end

  function RiskPriors{M}(ϵ::F, Ωs::F, Ωi::F, κ::F) where {
                         F <: Vector{UnivariateDistribution}, S <: SI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, nothing)
  end
end

function _indices(x::RiskPriors{M};
                  zeros::Bool=true) where {
                  S <: DiseaseStateSequence,
                  M <: ILM{S}}
  indices = [length(x.sparks)
             length(x.susceptibility)
             length(x.infectivity)
             length(x.transmissibility)]
  if S in [SEIR; SEI]
    push!(indices, length(x.latency))
  elseif zeros
    push!(indices, 0)
  end
  if S in [SEIR; SIR]
    push!(indices, length(x.removal))
  elseif zeros
    push!(indices, 0)
  end
  return cumsum(indices)
end

function Base.length(x::RiskPriors{M}) where M <: ILM
  return _indices(x)[end]
end

function Base.getindex(x::RiskPriors{M},
                       i::Int64) where M <: ILM
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end

function Base.setindex!(x::RiskPriors{M},
                        z::UnivariateDistribution,
                        i::Int64) where M <: ILM
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  getfield(x, riskfunc)[end - (indices[riskfunc] - i)] = z
  return x
end

function Base.show(io::IO, x::RiskPriors{M}) where M <: ILM
  return print(io, "$M risk function priors")
end
