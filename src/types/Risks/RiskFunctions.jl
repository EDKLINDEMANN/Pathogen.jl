struct RiskFunctions{M <: ILM}
  sparks::Union{Nothing, Function}
  susceptibility::Union{Nothing, Function}
  infectivity::Union{Nothing, Function}
  transmissibility::Union{Nothing, Function}
  latency::Union{Nothing, Function}
  removal::Union{Nothing, Function}

  function RiskFunctions{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F, Ωr::F) where {
                            F <: Function, S <: SEIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, Ωl)
  end

  function RiskFunctions{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωl::F) where {
                            F <: Function, S <: SEI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, Ωl, nothing)
  end

  function RiskFunctions{M}(ϵ::F, Ωs::F, Ωi::F, κ::F, Ωr::F) where {
                            F <: Function, S <: SIR, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, Ωr)
  end

  function RiskFunctions{M}(ϵ::F, Ωs::F, Ωi::F, κ::F) where {
                            F <: Function, S <: SI, M <: ILM{S}}
    return new(ϵ, Ωs, Ωi, κ, nothing, nothing)
  end
end

function Base.show(io::IO, x::RiskFunctions{M}) where M <: ILM
  return print(io, "$M risk functions")
end
