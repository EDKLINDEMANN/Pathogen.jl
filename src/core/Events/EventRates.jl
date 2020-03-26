struct EventRates{S <: DiseaseStateSequence}
  exposure::Union{Nothing, Vector{Float64}}
  infection::Union{Nothing, Vector{Float64}}
  removal::Union{Nothing, Vector{Float64}}
  individuals::Int64

  function EventRates{S}(n::Int64) where S <: SEIR
    return new{S}(fill(0.0, n), fill(0.0, n), fill(0.0, n), n)
  end

  function EventRates{S}(n::Int64) where S <: SEI
    return new{S}(fill(0.0, n), fill(0.0, n), nothing, n)
  end
  
  function EventRates{S}(n::Int64) where S <: SIR
    return new{S}(nothing, fill(0.0, n), fill(0.0, n), n)
  end  
  
  function EventRates{S}(n::Int64) where S <: SI
    return new{S}(nothing, fill(0.0, n), nothing, n)
  end
end

function Base.getindex(x::EventRates{S}, new_state::DiseaseState) where S <: DiseaseStateSequence
  if new_state == State_E && S <: Union{SEIR, SEI}
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R && S <: Union{SEIR, SIR}
    return x.removal
  else
    @error "Invalid indexing disease state"
  end
end

function Base.getindex(x::EventRates{S}, new_states::DiseaseStates) where S <: DiseaseStateSequence
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end

function Base.sum(x::EventRates{S}) where S <: DiseaseStateSequence
  return sum([sum(x[i]) for i in convert(DiseaseStates, S)[2:end]])
end