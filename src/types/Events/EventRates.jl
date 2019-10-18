struct EventRates{M <: ILM}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function EventRates{M}(n::Int64) where {S <: SEIR, M <: ILM{S}}
    return new{M}(fill(0.0, n),
                  fill(0.0, n),
                  fill(0.0, n),
                  n)
  end

  function EventRates{M}(n::Int64) where {S <: SEI, M <: ILM{S}}
    return new{M}(fill(0.0, n),
                  fill(0.0, n),
                  nothing,
                  n)
  end

  function EventRates{M}(n::Int64) where {S <: SIR, M <: ILM{S}}
    return new{M}(nothing,
                  fill(0.0, n),
                  fill(0.0, n),
                  n)
  end
  function EventRates{M}(n::Int64) where {S <: SI, M <: ILM{S}}
    return new{M}(nothing,
                  fill(0.0, n),
                  nothing,
                  n)
  end
end

function Base.getindex(x::EventRates{M}, new_state::DiseaseState) where M <: ILM
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    throw(BoundsError("Unrecognized indexing disease state"))
  end
end

function Base.getindex(x::EventRates{M}, new_states::Vector{DiseaseState}) where M <: ILM
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end

function Base.sum(x::EventRates{T}) where M <: ILM
  return sum([sum(x[i]) for i in _state_progressions[T][2:end]])
end
