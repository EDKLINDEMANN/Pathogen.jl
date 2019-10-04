mutable struct Events{T <: EpidemicModel}
  exposure::Vector{Float64}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64
  start_time::Float64

  function Events{T}(e::V, i::V, r::V; start_time::Float64=0.0) where {T <: SEIR, V <: Vector{Float64}, F <: Float64}
    if length(unique((length.([e; i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    return new{T}(e, i, r, length(i), start_time)
  end

  function Events{T}(n_ids::Int64; start_time::Float64=0.0) where T <: SEIR
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids), fill(NaN, n_ids), start_time=start_time)
  end

  function Events{T}(e::V, i::V; start_time::Float64=0.0) where {T <: SEI, V <: Vector{Float64}}
    if length(unique((length.([e; i])))) != 1
      @error "Length of event time vectors must be equal"
    end
    x = new{T}()
    x.exposure = e
    x.infection = i
    x.individuals = length(i)
    x.start_time = start_time
    return x
  end

  function Events{T}(n_ids::Int64; start_time::Float64=0.0) where T <: SEI
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids), start_time=start_time)
  end

  function Events{T}(i::V, r::V; start_time::Float64=0.0) where {T <: SIR, V <: Vector{Float64}}
    if length(unique((length.([i; r])))) != 1
      @error "Length of event time vectors must be equal"
    end
    x = new{T}()
    x.infection = i
    x.removal = r
    x.individuals = length(i)
    x.start_time = start_time
    return x
  end

  function Events{T}(n_ids::Int64; start_time::Float64=0.0) where T <: SIR
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids), start_time = start_time)
  end

  function Events{T}(i::V; start_time::Float64=0.0) where {T <: SI, V <: Vector{Float64}}
    x = new{T}()
    x.infection = i
    x.individuals = length(i)
    x.start_time = start_time
    return x
  end

  function Events{T}(n_ids::Int64; start_time::Float64=0.0) where T <: SI
    return Events{T}(fill(NaN, n_ids), start_time = start_time)
  end

  function Events{T}(a::Array{Float64,2}; start_time::Float64=0.0) where T <: EpidemicModel
    if size(a, 2) == 4
      return Events{T}(a[:,1], a[:,2], a[:,3], a[:,4], start_time = start_time)
    elseif size(a, 2) == 3
      return Events{T}(a[:,1], a[:,2], a[:,3], start_time = start_time)
    elseif size(a, 2) == 2
      return Events{T}(a[:,1], a[:,2], start_time = start_time)
    else
      @error "Invalid array size for construction of an $(Events{T}) object"
    end
  end

  function Events{T}(x::Vector{DiseaseState}; start_time::Float64=0.0) where T <: EpidemicModel
    events = Events{T}(length(x), start_time = start_time)
    for i = 1:length(x)
      if x[i] in _state_progressions[T]
        for j = _state_progressions[T][2:findfirst(Ref(x[i]) .== _state_progressions[T])]
          events[j][i] = -Inf
        end
      else
        @error "Invalid initial state for individual $i"
      end
    end
    return events
  end
end

function Base.show(io::IO, x::Events{T}) where T <: EpidemicModel
  return print(io, "$T model event times (n=$(x.individuals))")
end

function Base.getindex(x::Events{T}, new_state::DiseaseState) where T <: EpidemicModel
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    @error "Unrecognized indexing disease state"
  end
end

function Base.getindex(x::Events{T}, states::Vector{DiseaseState}) where T <: EpidemicModel
  y = x[states[1]]
  for i = 2:length(states)
    y = hcat(y, x[states[i]])
  end
  return y
end

function Base.convert(::Type{Array{Float64, 2}}, x::Events{T}) where T <: EpidemicModel
  return x[_state_progressions[T][2:end]]
end

function Base.convert(::Type{Vector{Float64}}, x::Events{T}) where T <: EpidemicModel
  return x[_state_progressions[T][2:end]][:]
end

function Base.convert(::Type{Array{Float64, 2}}, x::Array{Events{T}, 1}) where T <: EpidemicModel
  y = convert(Vector{Float64}, x[1])'
  for i = 2:length(x)
    y = vcat(y, convert(Vector{Float64}, x[i])')
  end
  return y
end

function Base.minimum(x::Events{T}) where T <: EpidemicModel
  y = convert(Array{Float64, 2}, x)
  return minimum(y[.!isnan.(y)])
end

function Base.maximum(x::Events{T}) where T <: EpidemicModel
  y = convert(Array{Float64, 2}, x)
  return maximum(y[.!isnan.(y)])
end

function Statistics.mean(x::Vector{Events{T}}) where T <: EpidemicModel
  return Events{T}(mean([convert(Array{Float64, 2}, i) for i in x]))
end
