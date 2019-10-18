struct RiskParameters{M <: ILM}
  sparks::Union{Nothing, Vector{Float64}}
  susceptibility::Union{Nothing, Vector{Float64}}
  infectivity::Union{Nothing, Vector{Float64}}
  transmissibility::Union{Nothing, Vector{Float64}}
  latency::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}

  function RiskParameters{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F, Ωr::F) where {
                            F <: Vector{Float64}, S <: SEIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, Ωl)
  end

  function RiskParameters{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F) where {
                            F <: Vector{Float64}, S <: SEI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, nothing)
  end

  function RiskParameters{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωr::F) where {
                            F <: Vector{Float64}, S <: SIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, Ωr)
  end

  function RiskParameters{M}(ϵ::F, Ωs::F, Ωi::F, κ::F) where {
                            F <: Vector{Float64}, S <: SI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, nothing)
  end
end

function _indices(x::RiskParameters{M};
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

function Base.length(x::RiskParameters{M}) where M <: ILM
  return _indices(x)[end]
end

function Base.getindex(x::RiskParameters{M},
                       i::Int64) where M <: ILM
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  return getfield(x, riskfunc)[end - (indices[riskfunc] - i)]
end

function Base.setindex!(x::RiskParameters{M},
                        z::Float64,
                        i::Int64) where M <: ILM
  indices = _indices(x, zeros = true)
  riskfunc = findfirst(i .<= indices)
  getfield(x, riskfunc)[end - (indices[riskfunc] - i)] = z
  return x
end

function Base.convert(::Type{Vector{Float64}},
                      x::RiskParameters{M}) where M <: ILM
  return [x[i] for i = 1:length(x)]
end

function Base.convert(::Type{Array{Float64, 2}},
                      x::Vector{RiskParameters{M}}) where M <: ILM
  return [x[i][j] for i = 1:length(x), j = 1:length(x[1])]
end

function _like(x::RiskParameters{M},
               v::Vector{Float64}) where {
               S <: DiseaseStateSequence,
               M <: ILM{S}}
  indices = _indices(x, zeros=false)
  if indices[end] != length(v)
    @error "Incompatiable parameter vector"
  end
  if S == SEIR
    return RiskParameters{M}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])],
                             v[(indices[5]+1):(indices[6])])
  elseif S in [SEI; SIR]
    return RiskParameters{M}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])],
                             v[(indices[4]+1):(indices[5])])
  elseif S == SI
    return RiskParameters{M}(v[1:(indices[1])],
                             v[(indices[1]+1):(indices[2])],
                             v[(indices[2]+1):(indices[3])],
                             v[(indices[3]+1):(indices[4])])
  end
end

function Base.show(io::IO, x::RiskParameters{M}) where M <: ILM
  return print(io, "$M risk function parameters")
end
